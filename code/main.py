from scipy.interpolate import griddata
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import basemap
import os.path
import util_fcns
from util_fcns import *
import pickle
import sys
import matplotlib as mpl

########################## INPUT and SETUP ###########################################

from da_config import *

vars_3d = ('XLONG','XLAT','HGT','T2','U10','V10','PBLH','LANDMASK','RAINNC','SNOWNC')
vars_4d = ('T','PH','PHB','P','PB','QRAIN')


if len(sys.argv) !=2:
	sys.exit('Please provide source NetCDF file as argument')

netcdf_name = sys.argv[1]
year, month, day, hr = netcdf_name[-19:-15], netcdf_name[-14:-12], netcdf_name[-11:-9], netcdf_name[-8:-6]
timestamp = netcdf_name[-19:-6]
plot_timestamp = '%s00UTC %s-%s-%s' %(hr,day,month,year)

#naming for saved plots
domain_plot = 'DEM_obs_12arcsec_%s.pdf' %timestamp
temp_plot = 'Downscaled_T2_%s.pdf' %timestamp
verif_plot = 'Downscaled_verif_T2_%s.pdf' %timestamp
corrected_verif_plot = 'Corrected_verif_T2_%s.pdf' %timestamp
corrected_temp_plot = 'Corrected_T2_%s.pdf' %timestamp
landmask_name = 'landmask_%s_%s_%s_%s.npy' %(lon[0],lon[1],lat[0],lat[1])
basemap_name = 'basemap_%s_%s_%s_%s_f.pickle'  %(lon[0],lon[1],lat[0],lat[1])
fig_subdir = os.path.join(fig_dir,year,month,day,hr)+'/'

########################### end of input ##########################################

#==================================SECTION 1: domain and data prep=============================
#get model data
wrf_data = os.path.join(netcdf_dir,year,month,day,netcdf_name)
fc_data = import_wrf_data(wrf_data,vars_3d,vars_4d)

#get with DEM data
elevation, lat_grid, lon_grid, bounds = import_dem_data(geotiff,lat,lon)

#configure basemaps
basemap_path = basemap_dir + basemap_name
if os.path.isfile(basemap_path):
	bm = pickle.load(open(basemap_path,'rb'))   # load here the above pickle 
	print('Domain basemap found at: %s' %basemap_path)
else:
	print('WARNING: no existing basemaps found: configuring a new basemap')
	bm = basemap.Basemap(llcrnrlon=lon[0], llcrnrlat=lat[0], urcrnrlon=lon[1], urcrnrlat=lat[1], resolution='f')
	pickle.dump(bm,open(basemap_path,'wb'),-1)  	# pickle the new map for later
	print('......... New basemap instance saved as: %s' %basemap_path)

#extract locations of observation stations
obs, obs_len = import_emwxnet_data(emx_name,os.path.join(emx_dir,year,month,day),bounds,hr)

#==================================SECTION 2:Prepare for Interpolation=============================

#load landmask. If unavailable generate new one (ATTN: SLOW ROUTINE!!!!)
landmask_path = landmask_dir + landmask_name
if os.path.isfile(landmask_path):
	dem_landmask = np.load(landmask_path)
	print('Domain landmask found at: %s' %landmask_path)
else:
	print('WARNING: No existing landmask found. Generating new landmask file: THIS MAY TAKE SEVERAL MINUTES')
	dem_landmask = make_landmask(bm,lon_grid,lat_grid)
	bm_fig = plot_domain(bm)
	bm.imshow(dem_landmask.T, origin='upper')
	bm.drawcoastlines()
	plt.savefig(fig_subdir + 'landmask_%s.pdf' %timestamp, format='pdf')
	plt.close()
	print('......... Domain landmask generated from basemaps saved as: %s' %landmask_path)
	np.save(landmask_path, dem_landmask)

print('Performing preliminary interpolation of model data')
fcx,fcy = bm(fc_data['XLONG'][:,:], fc_data['XLAT'][:,:])
fcPoints = zip(fcx.ravel(), fcy.ravel())
fcH, fcT, fcR = fc_data['HGT'][:,:], fc_data['T2'][:,:], fc_data['RAINNC'][:,:]
fcT_flat, fcH_flat, fcR_flat  = fcT.ravel(), fcH.ravel(), fcR.ravel()

#build KDTree and clean up obs
obs, demTree = prep_da(obs,lon_grid,lat_grid,elevation)
# np.savetxt('9h_delay.txt', obs['id'], fmt='%d')


#=============================SECTION 3: Temperature downscaling and DA=============================

interpT = split_interp(fc_data,lon_grid,lat_grid,dem_landmask)
interpH = griddata(fcPoints, fcH_flat, (lon_grid, lat_grid), method='cubic')

print('Calculating model lapse rate field')
fcGamma_flat = get_gamma(fc_data['T'],fc_data['P'],fc_data['PB'],fc_data['PH'],fc_data['PHB'],lvl)
interpGamma = griddata(fcPoints, fcGamma_flat, (lon_grid, lat_grid), method='cubic')

#do lapse rate adjustment 
print('Downscaling model T data')
demT = lapse_adjust(interpT,elevation,interpH,interpGamma)

#perform basic verification and produce scatter plots
print('Calculating T downscaling evaluation statistics')
verif_fig = verif_sets(fcPoints,demTree,obs,fcT,demT,'downscaling',plot_timestamp,dem_landmask,land_test=True,output='fig')
plt.savefig(fig_subdir+verif_plot, format='pdf')
print('......... Verification plots for downscaled T saved as: %s' %verif_plot)
plt.close()


print('Assimilating observational data for T')
#split observations into training and testing datasets (cross-evaluation)
obsTrain, obsTest = subset_obs(obs,verif_frac)

if run_MD:
	final_adjusted_T = da_md_land(obsTrain,elevation,lon_grid,lat_grid,dem_landmask,demT,interpGamma,params,bias_mode,dist_cutoff,bm)
	plot_tag = 'MD DA'
else:
	final_adjusted_T = da_roi(roi,elev_roi,obsTrain,elevation,demT,interpGamma,lon_grid,lat_grid,dem_landmask,bias_mode)
	plot_tag = 'Inverse distance DA'



print('Performing cross-evaluation of corrected T data')
verif_fig = verif_sets(fcPoints,demTree,obsTest,fcT,final_adjusted_T,plot_tag,plot_timestamp,dem_landmask,land_test=True,output='fig')
plt.savefig(fig_subdir+corrected_verif_plot, format='pdf')
print('......... Verification plots for corrected T datasets saved as: %s' %corrected_verif_plot)
plt.close()


#==================================SECTION 4: Precipitation Interpolaton================================

# interpR = griddata(fcPoints, fcR_flat, (lon_grid, lat_grid), method='cubic')

# bm_fig = plot_domain(bm)
# bm.imshow(elevation.T, cmap=plt.cm.gist_yarg, alpha = 0.7, origin='upper' )
# bm.drawcoastlines()
# plt.title('ITERPOLTED ACCUMULATED RAIN FIELD', fontweight='bold')
# bm.imshow(interpR.T, origin='upper', alpha = 0.7, cmap=plt.cm.ocean_r)
# cbar = plt.colorbar(label='MD weight',orientation='horizontal')
# cbar.solids.set_edgecolor('face')
# plt.show()


# plt.scatter(fcR_flat,fcH_flat)
# plt.title('PRECIP vs ELEVATION: %s' %plot_timestamp)
# plt.xlabel('Rainfall (mm)')
# plt.ylabel('Elevation (m)')
# plt.show()

# #collapse Qrain column and interpolate
# fcPLR = fc_data['QRAIN']
# fcPLR[fcPLR==0] = np.nan
# fcPLR = np.nanmean(fcPLR, axis = 0)
# fcPLR[np.isnan(fcPLR)] = 0
# fcPLR_flat = fcPLR.ravel()
# interpPLR = griddata(fcPoints, fcPLR_flat, (lon_grid, lat_grid), method='cubic')

# bm_fig = plot_domain(bm)
# bm.imshow(elevation.T, cmap=plt.cm.gist_yarg, alpha = 0.7, origin='upper' )
# bm.drawcoastlines()
# # plt.title('Station Locations and MD Weights', fontweight='bold')
# bm.imshow(interpPLR.T, origin='upper', alpha = 0.7, cmap=plt.cm.ocean_r)
# cbar = plt.colorbar(label='MD weight',orientation='horizontal')
# cbar.solids.set_edgecolor('face')
# plt.show()

# #plot vertical slice of Qr and elevation
# fig, ax1 = plt.subplots()
# ax1.contourf(fcPLR[:,360,:])
# ax2 = ax1.twinx()
# ax1.set_ylabel('Qrain', color='b')
# ax2.plot(fc_data['HGT'][360,:],'r')
# ax2.set_ylabel('elevation (m)', color='r')
# plt.show()

# final_adjusted_R = da_md(obsTrain,elevation,lon_grid,lat_grid,interpR,params,bias_mode,dist_cutoff)
# verif_fig = verif_sets(fcPoints,demTree,obsTest,fcR,final_adjusted_R,plot_tag,plot_timestamp,dem_landmask,land_test=True,output='fig')
# plt.show()

# bm_fig = plot_domain(bm)
# bm.drawcoastlines()
# # plt.title('Station Locations and MD Weights', fontweight='bold')
# im = bm.imshow(final_adjusted_R.T, origin='upper', alpha = 0.7, cmap=plt.cm.ocean_r)
# im.set_clim([0,10])
# cbar = plt.colorbar(label='MD weight',orientation='horizontal')
# cbar.solids.set_edgecolor('face')
# plt.show()


#===============================SECTION 5: Saving Data and Generating plots=============================

#setup colormaps
load_cmap = np.load('./npy/Tcmap.npy')
hot = mpl.colors.ListedColormap(load_cmap/255)
cool = plt.cm.cubehelix
cool_vals = [cool(i) for i in range(cool.N)]
hot_vals = [hot(i) for i in range(hot.N)]
comb_vals = cool_vals + hot_vals
cmT = mpl.colors.ListedColormap(comb_vals)


print('Plotting useful data')
#plot downscaled temperature
bm_fig = plot_domain(bm)
im = bm.imshow(demT.T-273, origin='upper', cmap=cmT)
im.set_clim(T_range)
cbar = plt.colorbar(label='temperature[C]',orientation='horizontal')
cbar.solids.set_edgecolor('face')
plt.title("DOWNSCALED 2m TEMPERATURE | %s" %plot_timestamp, fontweight='bold')
plt.savefig(fig_subdir+temp_plot, format='pdf')
print('......... Downscaled 2m temperature fields saved as: %s' %temp_plot)
plt.close()

#plot DA-corrected temperature
bm_fig = plot_domain(bm)
im = bm.imshow(final_adjusted_T.T-273, origin='upper', cmap=cmT)
im.set_clim(T_range)
cbar = plt.colorbar(label='temperature[C]',orientation='horizontal')
cbar.solids.set_edgecolor('face')
scat = bm.scatter(obsTrain['x'],obsTrain['y'],linewidth='0.5', s=17,marker='o', c=np.array(obsTrain['t']-273), cmap=cmT)
scat.set_clim(T_range)
plt.title("HIGH-RESOLUTION TEMPERATURE ANALYSIS (2m) | %s | %s" %(plot_timestamp,plot_tag))
plt.savefig(fig_subdir+corrected_temp_plot, format='pdf')
print('......... Corrected T at 2m temperature fields saved as: %s' %corrected_temp_plot)
plt.close()

#plot observation stations overlayed on topo
plt.figure(figsize=(30, 15))
bm.drawcoastlines(linewidth=0.3,color='grey')
plt.title('DEM MODEL WITH SELECT EMWXNET STATIONS | %s' %plot_timestamp, fontweight='bold')
bm.imshow(elevation.T, cmap=plt.cm.gist_earth, origin='upper')
cbar = plt.colorbar(label='elevation [m]', orientation='horizontal')
cbar.solids.set_edgecolor('face')
bm.scatter(obs['x'],obs['y'],s=6,marker='o', color ='r')
plt.savefig(fig_subdir+domain_plot, format='pdf')
plt.close()
print('......... Domain with DEM model and observations saved as: %s' %domain_plot)

#plot lapse rate
plt.figure(figsize=(30, 15))
plt.title('INTERPOLATED MODEL LAPSE RATE | %s' %plot_timestamp, fontweight='bold')
bm.drawcoastlines(linewidth=0.3,color='grey')
im = bm.imshow(interpGamma.T, origin='upper')
im.set_clim([-8.5,-4.5])
cbar = plt.colorbar(label='lapse rate [deg C/km]',orientation='horizontal')
cbar.solids.set_edgecolor('face')
plt.savefig(fig_subdir+'lapse_rate_%s.pdf'%timestamp, format='pdf')
plt.close()
print('......... Interpolated model lapse rate saved as: %s' %('lapse_rate_%s.pdf' %timestamp))

#plot observations elevation histogram
plt.figure(figsize=(20, 10))
plt.title('HISTROGRAM OF EMWXNET STATION ELEVATIONS | %s' %plot_timestamp, fontweight='bold', fontsize=10)
plt.hist(obs['h'])
plt.xlabel('elevation (m)', fontsize=10)
plt.ylabel('station count', fontsize=10)
plt.savefig(fig_subdir+'obs_histogram_%s.pdf'%timestamp, format='pdf')
plt.close()
print('......... Observations histogram saved as: %s' %('obs_histogram_%s.pdf' %timestamp))

#plot raw temperature field
plt.figure(figsize=(30, 15))
plt.title('Raw model temperature field (2m): %s' %plot_timestamp, fontweight='bold')
bm.drawcoastlines(linewidth=0.3,color='grey')
im = bm.contourf(fcx,fcy, fcT-273, levels=np.arange(T_range[0],T_range[1],0.1), cmap=cmT)
im.set_clim(T_range)
cbar = plt.colorbar(label='temperature [C]',orientation='horizontal')
cbar.solids.set_edgecolor('face')
plt.savefig(fig_subdir+'temp_raw_%s.pdf' %timestamp, format='pdf')
plt.close()
print('......... Raw temperature fields saved as: %s' %('temp_raw_%s.pdf' %timestamp))

print('Run COMPLETE.')

# NOTES FOR FUTURE OPTIMIZATION
# -lapse_adjust - avoid point by point (work with arrays instead)
