from scipy.interpolate import griddata
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import basemap
import os.path
from util_fcns import *
import pickle
import sys
import matplotlib as mpl

########################## INPUT and SETUP ###########################################

from da_config import *

vars_3d = ('XLONG','XLAT','HGT','T2','U10','V10','PBLH','LANDMASK')
vars_4d = ('T','PH','PHB','P','PB')


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

#==================================SECTION 2: Temperature Interpolation=============================

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
fcH, fcT = fc_data['HGT'][:,:], fc_data['T2'][:,:]
fcT_flat, fcH_flat  = fcT.ravel(), fcH.ravel()

#build KDTree and clean up obs
obs, demTree = prep_da(obs,lon_grid,lat_grid,elevation)
# np.savetxt('9h_delay.txt', obs['id'], fmt='%d')

interpT = split_interp(fc_data,lon_grid,lat_grid,dem_landmask)
interpH = griddata(fcPoints, fcH_flat, (lon_grid, lat_grid), method='cubic')

print('Calculating model lapse rate field')
fcGamma_flat = get_gamma(fc_data['T'],fc_data['P'],fc_data['PB'],fc_data['PH'],fc_data['PHB'],lvl)
interpGamma = griddata(fcPoints, fcGamma_flat, (lon_grid, lat_grid), method='cubic')

#do lapse rate adjustment 
print('Downscaling model data')
demT = lapse_adjust(interpT,elevation,interpH,interpGamma)

#perform basic verification and produce scatter plots
print('Calculating downscaling evaluation statistics')
verif_fig = verif_sets(fcPoints,demTree,obs,fcT,demT,'downscaling',plot_timestamp,dem_landmask,land_test=True,output='fig')
plt.savefig(fig_subdir+verif_plot, format='pdf')
print('......... Verification plots saved as: %s' %verif_plot)
plt.close()

#==================================SECTION 3: Temperature Data Assimilation=============================

print('Assimilating observational data')

#split observations into training and testing datasets (cross-evaluation)
obsTrain, obsTest = subset_obs(obs,verif_frac)

if run_MD:
	final_adjusted_T = da_md(obsTrain,elevation,lon_grid,lat_grid,dem_landmask,demT,interpGamma,params,bias_mode,dist_cutoff,bm)
	plot_tag = 'MD DA'
else:
	final_adjusted_T = da_roi(roi,elev_roi,obsTrain,elevation,demT,interpGamma,lon_grid,lat_grid,dem_landmask,bias_mode)
	plot_tag = 'Inverse distance DA'



print('Performing cross-evaluation of corrected T data')
verif_fig = verif_sets(fcPoints,demTree,obsTest,fcT,final_adjusted_T,plot_tag,plot_timestamp,dem_landmask,land_test=True,output='fig')
plt.savefig(fig_subdir+corrected_verif_plot, format='pdf')
print('......... Verification plots for corrected datasets saved as: %s' %corrected_verif_plot)
plt.close()


#===============================SECTION 4: Saving Data and Generating plots=============================

cool = plt.cm.cubehelix
# hot = plt.cm.CMRmap_r
load_cmap = np.load('./npy/Tcmap.npy')
cmT = mpl.colors.ListedColormap(load_cmap/255)
hot = cmT

cool_vals = [cool(i) for i in range(cool.N)]
hot_vals = [hot(i) for i in range(hot.N)]

comb_vals = cool_vals + hot_vals

# random hue with constant sat and value
cmT = mpl.colors.ListedColormap(comb_vals)
# load_cmap = np.load('./npy/Tcmap.npy')
# cmT = mpl.colors.ListedColormap(load_cmap/255)

print('Plotting useful data')
#plot downscaled temperature
bm_fig = plot_domain(bm)
# im = bm.imshow(demT.T-273, origin='upper', cmap=plt.cm.gist_ncar, alpha=.7)
# im = bm.imshow(demT.T-273, origin='upper', cmap=plt.cm.Spectral_r)
im = bm.imshow(demT.T-273, origin='upper', cmap=cmT)
im.set_clim(T_range)
cbar = plt.colorbar(label='temperature[C]',orientation='horizontal')
cbar.solids.set_edgecolor('face')
# bm.contour(lon_grid,lat_grid,demT-273, colors='w', alpha=0.5, linewideth=0.05)
plt.title("DOWNSCALED 2m TEMPERATURE | %s" %plot_timestamp, fontweight='bold')
plt.savefig(fig_subdir+temp_plot, format='pdf')
print('......... Downscaled 2m temperature fields saved as: %s' %temp_plot)
plt.close()

#plot DA-corrected temperature
bm_fig = plot_domain(bm)
# im = bm.imshow(final_adjusted_T.T-273, origin='upper', cmap=plt.cm.gist_ncar, alpha=0.7)
# im = bm.imshow(final_adjusted_T.T-273, origin='upper', cmap=plt.cm.Spectral_r)
im = bm.imshow(final_adjusted_T.T-273, origin='upper', cmap=cmT)
im.set_clim(T_range)
cbar = plt.colorbar(label='temperature[C]',orientation='horizontal')
cbar.solids.set_edgecolor('face')
bm.scatter(obsTrain['x'],obsTrain['y'],s=6,marker='o', color ='k')
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
print('......... Interpolated model lapse rate saved as: %s' %(fig_dir +'lapse_rate_%s.pdf' %timestamp))

#plot observations elevation histogram
plt.figure(figsize=(20, 10))
plt.title('HISTROGRAM OF EMWXNET STATION ELEVATIONS | %s' %plot_timestamp, fontweight='bold', fontsize=10)
plt.hist(obs['h'])
plt.xlabel('elevation (m)', fontsize=10)
plt.ylabel('station count', fontsize=10)
plt.savefig(fig_subdir+'obs_histogram_%s.pdf'%timestamp, format='pdf')
plt.close()
print('......... Observations histogram saved as: %s' %(fig_dir +'obs_histogram_%s.pdf' %timestamp))

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
print('......... Raw temperature fields saved as: %s' %(fig_dir +'temp_raw_%s.pdf' %timestamp))

print('Run COMPLETE.')




'''
TESTED LAPSE RATE ROUTINES:
-based on surface forecast forecast points 
-based on obs surface points 
-based on model level data  - best result, most robust


OPTIMIZATION POSSIBILITIES:
-avoid rebuilding kd-tree for MD
-lapse_adjust - avoid point by point (work with arrays instead)
'''
