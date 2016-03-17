#!/usr/bin/env python
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
from matplotlib.colors import LogNorm

########################## INPUT and SETUP ###########################################

from da_config import *

vars_3d = ('XLONG','XLAT','HGT','T2','PBLH','U10','V10','LANDMASK','RAINNC','SNOWNC')
vars_4d = ('T','PH','PHB','P','PB','QRAIN')

if len(sys.argv) !=2:
	sys.exit('Please provide source NetCDF file as argument')

netcdf_name = sys.argv[1]
year, month, day, hr = netcdf_name[-19:-15], netcdf_name[-14:-12], netcdf_name[-11:-9], netcdf_name[-8:-6]
timestamp = netcdf_name[-19:-6]
plot_timestamp = '%s00UTC %s-%s-%s' %(hr,day,month,year)

#naming for saved plots
domain_plot = 'DEM_obs_%s.pdf' %timestamp
ds_plot = 'downscaled_%s.pdf' %timestamp
hs_plot = 'highres_%s.pdf' %timestamp
ds_verif_plot = 'downscaled_verif_%s.pdf' %timestamp
hs_verif_plot = 'highres_verif_%s.pdf' %timestamp
landmask_name = 'landmask_%s_%s_%s_%s.npy' %(lon[0],lon[1],lat[0],lat[1])
basemap_name = 'basemap_%s_%s_%s_%s_f.pickle'  %(lon[0],lon[1],lat[0],lat[1])
fig_subdir = os.path.join(fig_dir,year,month,day,hr)+'/'

########################### end of input ##########################################

#==================================SECTION 1: domain and data prep=============================
#get model data
wrf_data = os.path.join(netcdf_dir,year,month,day,netcdf_name)
fc_data = import_wrf_data(wrf_data,vars_3d,vars_4d)

#get with DEM data
elevation, lat_grid, lon_grid, bounds = import_dem_data(elev_geotiff,lat,lon)

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
obsT, obsT_len, obsW, obsW_len, obsR, obsR_len = import_emwxnet_data(emx_name,os.path.join(emx_dir,year,month,day),bounds,hr)

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

print('Preparing model data')
fcx,fcy = bm(fc_data['XLONG'][:,:], fc_data['XLAT'][:,:])
fcPoints = zip(fcx.ravel(), fcy.ravel())
fcH, fcT, fcR, fcU, fcV = fc_data['HGT'][:,:], fc_data['T2'][:,:], fc_data['RAINNC'][:,:], fc_data['U10'][:,:], fc_data['V10'][:,:]
fcH, fcT, fcR, fcU, fcV = fc_data['HGT'][:,:], fc_data['T2'][:,:], fc_data['RAINNC'][:,:], fc_data['U10'][:,:], fc_data['V10'][:,:]
fcT_flat, fcH_flat, fcR_flat, fcU_flat, fcV_flat  = fcT.ravel(), fcH.ravel(), fcR.ravel(), fcU.ravel(), fcV.ravel()

#build KDTree and clean up obs
print('Testing observations metadata: Temperature')
obsT, demTreeT = prep_da(obsT,lon_grid,lat_grid,elevation)
print('Testing observations metadata: Wind')
obsW, demTreeW = prep_da(obsW,lon_grid,lat_grid,elevation)
print('Testing observations metadata: Rain')
obsR, demTreeR = prep_da(obsR,lon_grid,lat_grid,elevation)

# =============================SECTION 3: Temperature downscaling and DA=============================

interpT = split_interp(fc_data,lon_grid,lat_grid,dem_landmask,var='T2')
interpH = griddata(fcPoints, fcH_flat, (lon_grid, lat_grid), method='cubic')

print('Calculating model temperature lapse rate field')
fcGamma_flat = get_gamma(fc_data['T'],fc_data['P'],fc_data['PB'],fc_data['PH'],fc_data['PHB'],lvl)
interpGamma = griddata(fcPoints, fcGamma_flat, (lon_grid, lat_grid), method='cubic')

#do lapse rate adjustment 
print('Downscaling model T data')
demT = lapse_adjust(interpT,elevation,interpH,interpGamma)

#perform basic verification and produce scatter plots
print('Calculating T downscaling evaluation statistics')
verif_fig = verif_sets(fcPoints,demTreeT,obsT,fcT,demT,'downscaling',plot_timestamp,dem_landmask,land_test=True,var='t')
plt.savefig(fig_subdir+'T2_'+ds_verif_plot, format='pdf')
print('......... Verification plots for downscaled T saved as: %s' %('T2_'+ds_verif_plot))
plt.close()

print('Assimilating observational data for T')
#split observations into training and testing datasets (cross-evaluation)
obsTrainT, obsTestT = subset_obs(obsT,verif_frac)

if temp_run_MD:
	final_adjusted_T = da_md(obsTrainT,elevation,lon_grid,lat_grid,dem_landmask,demT,params,temp_bias_mode,dist_cutoff,bm,interpGamma,land_mode=1,var='t')
	plot_tag = 'MD DA'
else:
	final_adjusted_T = da_roi(temp_roi,temp_elev_roi,obsTrainT,elevation,demT,lon_grid,lat_grid,dem_landmask,temp_bias_mode,interpGamma,var='t')
	plot_tag = 'Inverse distance DA'

print('Performing cross-evaluation of corrected T data')
verif_fig = verif_sets(fcPoints,demTreeT,obsTestT,fcT,final_adjusted_T,plot_tag,plot_timestamp,dem_landmask,land_test=True,var='t')
plt.savefig(fig_subdir+'T2_'+hs_verif_plot, format='pdf')
print('......... Verification plots for corrected T datasets saved as: %s' %('T2_'+hs_verif_plot))
plt.close()

# #================================SECTION 4: Wind Interpolaton and DA================================

# print('Interpolating model wind data')
# # demU = split_interp(fc_data,lon_grid,lat_grid,dem_landmask,var='U10')
# # demV = split_interp(fc_data,lon_grid,lat_grid,dem_landmask,var='V10')
# interpU = split_interp(fc_data,lon_grid,lat_grid,dem_landmask,var='U10')
# interpV = split_interp(fc_data,lon_grid,lat_grid,dem_landmask,var='V10')

# #import aspect data and check that it corresponds to dem
# aspect, aspect_lat_grid, aspect_lon_grid, aspect_bounds = import_dem_data(aspect_geotiff,lat,lon)
# slope, slope_lat_grid, slope_lon_grid, slope_bounds = import_dem_data(slope_geotiff,lat,lon)
# if not ((aspect_lon_grid==lon_grid).all() and (aspect_lat_grid==lat_grid).all()):
# 	sys.exit('ERROR: Mismatch between DEM and ASPECT geotiffs - aborting data assimilation ')
# if not ((slope_lon_grid==lon_grid).all() and (slope_lat_grid==lat_grid).all()):
# 	sys.exit('ERROR: Mismatch between DEM and SLOPE geotiffs - aborting data assimilation ')


# #perform 2D diagnostic downscaling
# demU, demV = np.empty_like(interpU)* np.nan, np.empty_like(interpU)* np.nan
# slope = slope/100.
# interpSp = np.hypot(interpU, interpV)
# r2d = 180./np.pi
# interpDir = np.arctan2(interpU, interpV) * r2d
# interpDir = (interpDir + 360.) % 360;
# plt.imshow(interpDir.T)
# plt.colorbar()
# plt.show()


# aspect[aspect<0] = np.nan
# aspect_pve = (aspect + 90.) %360
# aspect_nve = (aspect - 90.) %360
# diff_pve =  180 - abs(abs(interpDir - aspect_pve) - 180)
# diff_nve = 180 - abs(abs(interpDir - aspect_nve) - 180)

# projection_plane = np.empty_like(interpDir) * np.nan
# projection_plane[diff_pve<diff_nve] = diff_pve[diff_pve<diff_nve]
# projection_plane[diff_nve<diff_pve] = diff_nve[diff_nve<diff_pve]
# projection_plane[np.isnan(projection_plane)] = interpDir[np.isnan(projection_plane)]

# d2r = np.pi/180.
# demU = interpSp * np.sin(d2r*projection_plane) * slope + interpSp * np.sin(d2r*interpDir) * (1-slope)
# demV = interpSp * np.cos(d2r*projection_plane) * slope + interpSp * np.cos(d2r*interpDir) * (1-slope)

# # #####CONSIDERATION::::: difference between interpolated height and dem height - to weigh adjustment


# #plot vector wind field
# bm_fig = plot_domain(bm)
# M = np.hypot(demU, demV)
# bm.imshow(M.T, origin='upper', cmap=plt.cm.YlOrBr)
# im = bm.quiver(lon_grid[::30,::30],lat_grid[::30,::30], demU[::30,::30], demV[::30,::30], latlon=True, width=0.001)
# im.axes.set_aspect(1.7)
# cbar = plt.colorbar(label='wind speed [m/s]',orientation='horizontal')
# cbar.solids.set_edgecolor('face')
# plt.show()
# # plt.title("HIGH-RESOLUTION ANALYSIS VECTOR WIND FiELD | %s | %s" %(plot_timestamp,plot_tag))
# # plt.savefig(fig_subdir+'wind_'+hs_plot, format='pdf')
# # print('......... Corrected vector wind fields saved as: %s' %('wind_'+hs_plot))
# # plt.close()



# bm_fig = plot_domain(bm) 
# M = np.hypot(interpU, interpV)
# bm.imshow(M.T, origin='upper', cmap=plt.cm.YlOrBr)
# im = bm.quiver(lon_grid[::30,::30],lat_grid[::30,::30], interpU[::30,::30], interpV[::30,::30], latlon=True, width=0.001)
# im2 = bm.quiver(obsTrainW['x'],obsTrainW['y'],obsTrainW['u'],obsTrainW['v'],latlon=True, width=0.001, color='r')
# im.axes.set_aspect(1.7)
# cbar = plt.colorbar(label='wind speed [m/s]',orientation='horizontal')
# cbar.solids.set_edgecolor('face')
# plt.show()

# #perform basic verification and produce scatter plots
# print('Calculating wind interpolation evaluation statistics')
# verif_fig = verif_sets(fcPoints,demTreeW,obsW,fcU,demU,'interpolation',plot_timestamp,dem_landmask,land_test=False,var='u')
# plt.savefig(fig_subdir+'U10_'+ds_verif_plot, format='pdf')
# print('......... Verification plots for interpolated U saved as: %s' %('U10_'+ds_verif_plot))
# plt.close()
# verif_fig = verif_sets(fcPoints,demTreeW,obsW,fcV,demV,'interpolation',plot_timestamp,dem_landmask,land_test=False,var='v')
# plt.savefig(fig_subdir+'V10_'+ds_verif_plot, format='pdf')
# print('......... Verification plots for interpolated V saved as: %s' %('V10_'+ds_verif_plot))
# plt.close()

# print('Assimilating observational data for wind')
# #split observations into training and testing datasets (cross-evaluation)
# obsTrainW, obsTestW = subset_obs(obsW,verif_frac)

# if wind_run_MD:
# 	final_adjusted_U = da_md(obsTrainW,elevation,lon_grid,lat_grid,dem_landmask,demU,params,wind_bias_mode,dist_cutoff,bm,land_mode=0,var='u')
# 	final_adjusted_V = da_md(obsTrainW,elevation,lon_grid,lat_grid,dem_landmask,demV,params,wind_bias_mode,dist_cutoff,bm,land_mode=0,var='v')
# 	plot_tag = 'MD DA'
# else:
# 	final_adjusted_U = da_roi(wind_roi,wind_elev_roi,obsTrainW,elevation,demU,lon_grid,lat_grid,dem_landmask,wind_bias_mode,var='u')
# 	final_adjusted_V = da_roi(wind_roi,wind_elev_roi,obsTrainW,elevation,demV,lon_grid,lat_grid,dem_landmask,wind_bias_mode,var='v')
# 	plot_tag = 'Inverse distance DA'

# print('Performing cross-evaluation of corrected wind data')
# verif_fig = verif_sets(fcPoints,demTreeW,obsTestW,fcU,final_adjusted_U,plot_tag,plot_timestamp,dem_landmask,land_test=True,var='u')
# plt.savefig(fig_subdir+'U10_'+hs_verif_plot, format='pdf')
# print('......... Verification plots for corrected U datasets saved as: %s' %('U10_'+hs_verif_plot))
# plt.close()
# verif_fig = verif_sets(fcPoints,demTreeW,obsTestW,fcV,final_adjusted_V,plot_tag,plot_timestamp,dem_landmask,land_test=True,var='v')
# plt.savefig(fig_subdir+'V10_'+hs_verif_plot, format='pdf')
# print('......... Verification plots for corrected U datasets saved as: %s' %('V10_'+hs_verif_plot))
# plt.close()

# ==============================SECTION 5: Precipitation Interpolaton================================

demR = griddata(fcPoints, fcR_flat, (lon_grid, lat_grid), method='cubic')
obsTrainR, obsTestR = subset_obs(obsR,verif_frac)

if precip_run_PRISM:
	# final_adjusted_R = da_prism(obsTrainR,elevation,lon_grid,lat_grid,dem_landmask,demR,params,precip_bias_mode,dist_cutoff,bm,land_mode=0,var='r')
	final_adjusted_R = da_prism(obsTrainR,elevation,lon_grid,lat_grid,demR,int(month),lat,lon,rain_roi,bm,precip_bias_mode)
	plot_tag = 'PRISM DA'
else:
	final_adjusted_R = da_roi(rain_roi,precip_elev_roi,obsTrainR,elevation,demR,lon_grid,lat_grid,dem_landmask,precip_bias_mode,var='r')
	plot_tag = 'Inverse distance DA'

# final_adjusted_R = da_md(obsTrain,elevation,lon_grid,lat_grid,interpR,params,bias_mode,dist_cutoff)
print('Performing cross-evaluation of corrected wind data')
verif_fig = verif_sets(fcPoints,demTreeR,obsTestR,fcR,final_adjusted_R,plot_tag,plot_timestamp,dem_landmask,land_test=True,var='r')
plt.savefig(fig_subdir+'rain_'+hs_verif_plot, format='pdf')
print('......... Verification plots for corrected rainfall datasets saved as: %s' %('rain_'+hs_verif_plot))
plt.close()

# ===============================SECTION 5: Saving Data and Generating plots=============================

#-------------------plotting cosmetics-----------------
#generate colormaps
load_cmapT = np.load('./npy/Tcmap.npy')
hot = mpl.colors.ListedColormap(load_cmapT/255)
cool = plt.cm.cubehelix
cool_vals = [cool(i) for i in range(cool.N)]
hot_vals = [hot(i) for i in range(hot.N)]
comb_vals = cool_vals + hot_vals
cmT = mpl.colors.ListedColormap(comb_vals)

load_cmapR = np.load('./npy/Rcmap.npy')
cmR = mpl.colors.ListedColormap(load_cmapR/255)


#configure precip ranges for log plotting
if R_range[0]<=0:
	Rmin = 0.0001
else:
	Rmin = R_range[0]
Rmax = R_range[1]
final_adjusted_R_plot = np.copy(final_adjusted_R)
final_adjusted_R_plot[final_adjusted_R_plot<Rmin] = Rmin
rticks = range(10) + range(10,Rmax,10)


print('Plotting useful data')
#plot observation stations overlayed on topo
plt.figure(figsize=(15, 13))
bm.drawcoastlines(linewidth=0.3,color='grey')
plt.title('DEM MODEL WITH SELECT EMWXNET STATIONS | %s' %plot_timestamp, fontweight='bold')
im = bm.imshow(elevation.T, cmap=plt.cm.gist_earth, origin='upper')
im.axes.set_aspect(1.7)
cbar = plt.colorbar(label='elevation [m]', orientation='horizontal')
cbar.solids.set_edgecolor('face')
scat = bm.scatter(obsT['x'],obsT['y'],latlon=True, s=6,marker='o', color ='r')
plt.savefig(fig_subdir+domain_plot, format='pdf')
plt.close()
print('......... Domain with DEM model and observations saved as: %s' %domain_plot)

#plot downscaled temperature
bm_fig = plot_domain(bm)
im = bm.imshow(demT.T-273, origin='upper', cmap=cmT)
im.set_clim(T_range)
im.axes.set_aspect(1.7)
cbar = plt.colorbar(label='temperature[C]',orientation='horizontal')
cbar.solids.set_edgecolor('face')
plt.title("DOWNSCALED 2M TEMPERATURE | %s" %plot_timestamp, fontweight='bold')
plt.savefig(fig_subdir+'T2_'+ds_plot, format='pdf')
print('......... Downscaled 2m temperature fields saved as: %s' %('T2_'+ds_plot))
plt.close()

#plot DA-corrected temperature
bm_fig = plot_domain(bm)
im = bm.imshow(final_adjusted_T.T-273, origin='upper', cmap=cmT)
im.set_clim(T_range)
im.axes.set_aspect(1.7)
cbar = plt.colorbar(label='temperature[C]',orientation='horizontal')
cbar.solids.set_edgecolor('face')
scat = bm.scatter(obsTrainT['x'],obsTrainT['y'],linewidth='0.5', s=17,marker='o', c=np.array(obsTrainT['t']-273), cmap=cmT)
scat.set_clim(T_range)
<<<<<<< HEAD
plt.title("HIGH-RESOLUTION TEMPERATURE ANALYSIS (2M) | %s | %s" %(plot_timestamp,plot_tag))
plt.savefig(fig_subdir+'T2_'+hs_plot, format='pdf')
=======
plt.title("HIGH-RESOLUTION TEMPERATURE ANALYSIS (2M) | %s | %s" %(plot_timestamp,plot_tag_T))
bm_fig.tight_layout()
plt.savefig(fig_subdir+'T2_'+hs_plot, format='pdf', bbox_inches='tight')
>>>>>>> parent of 6f03a37... clean up main.py
print('......... Corrected T at 2m fields saved as: %s' %('T2_'+hs_plot))
plt.close()

#plot lapse rate
plt.figure(figsize=(15, 13))
plt.title('INTERPOLATED MODEL LAPSE RATE | %s' %plot_timestamp, fontweight='bold')
bm.drawcoastlines(linewidth=0.3,color='grey')
im = bm.imshow(interpGamma.T, origin='upper')
im.set_clim([-8.5,-4.5])
im.axes.set_aspect(1.7)
cbar = plt.colorbar(label='lapse rate [deg C/km]',orientation='horizontal')
cbar.solids.set_edgecolor('face')
plt.savefig(fig_subdir+'lapse_rate_%s.pdf'%timestamp, format='pdf')
plt.close()
print('......... Interpolated model lapse rate saved as: %s' %('lapse_rate_%s.pdf' %timestamp))

#plot observations elevation histogram
plt.figure(figsize=(20, 10))
plt.title('HISTROGRAM OF EMWXNET STATION ELEVATIONS | %s' %plot_timestamp, fontweight='bold', fontsize=10)
plt.hist(obsT['h'])
plt.xlabel('elevation (m)', fontsize=10)
plt.ylabel('station count', fontsize=10)
plt.savefig(fig_subdir+'obs_histogram_%s.pdf'%timestamp, format='pdf')
plt.close()
print('......... Observations histogram saved as: %s' %('obs_histogram_%s.pdf' %timestamp))


#plot raw temperature field
plt.figure(figsize=(15, 13))
plt.title('RAW MODEL TEMPERATURE FIELD (2M): %s' %plot_timestamp, fontweight='bold')
bm.drawcoastlines(linewidth=0.3,color='grey')
im = bm.contourf(fcx,fcy, fcT-273, levels=np.arange(T_range[0],T_range[1],0.1), cmap=cmT)
im.set_clim(T_range)
# im.axes.set_aspect(1.7)
cbar = plt.colorbar(label='temperature [C]',orientation='horizontal')
cbar.solids.set_edgecolor('face')
plt.savefig(fig_subdir+'temp_raw_%s.pdf' %timestamp, format='pdf')
plt.close()
print('......... Raw temperature fields saved as: %s' %('temp_raw_%s.pdf' %timestamp))


# #plot DA-corrected U wind component
# bm_fig = plot_domain(bm)	
# im = bm.imshow(final_adjusted_U.T, origin='upper', cmap=cmT)
# im.set_clim(W_range)
# im.axes.set_aspect(1.7)
# cbar = plt.colorbar(label='u- wind component',orientation='horizontal')
# cbar.solids.set_edgecolor('face')
# scat = bm.scatter(obsTrainW['x'],obsTrainW['y'],linewidth='0.5', s=17,marker='o', c=np.array(obsTrainW['u']), cmap=cmT)
# scat.set_clim(W_range)
# plt.title("HIGH-RESOLUTION U WIND COMPONENT ANALYSIS (10M) | %s | %s" %(plot_timestamp,plot_tag))
# plt.savefig(fig_subdir+'U10_'+hs_plot, format='pdf')
# print('......... Corrected U at 10M fields saved as: %s' %('U10_'+hs_plot))
# plt.close()


# #plot DA-corrected V wind component
# bm_fig = plot_domain(bm)
# im = bm.imshow(final_adjusted_V.T, origin='upper', cmap=cmT)
# im.set_clim(W_range)
# im.axes.set_aspect(1.7)
# cbar = plt.colorbar(label='v- wind component',orientation='horizontal')
# cbar.solids.set_edgecolor('face')
# scat = bm.scatter(obsTrainW['x'],obsTrainW['y'],linewidth='0.5', s=17,marker='o', c=np.array(obsTrainW['v']), cmap=cmT)
# scat.set_clim(W_range)
# plt.title("HIGH-RESOLUTION V WIND COMPONENT ANALYSIS (10M) | %s | %s" %(plot_timestamp,plot_tag))
# plt.savefig(fig_subdir+'V10_'+hs_plot, format='pdf')
# print('......... Corrected V at 10M fields saved as: %s' %('V10_'+hs_plot))
# plt.close()

# #plot vector wind field
# bm_fig = plot_domain(bm)
# M = np.hypot(final_adjusted_U, final_adjusted_V)
# bm.imshow(M.T, origin='upper', cmap=plt.cm.YlOrBr)
# im = bm.quiver(lon_grid[::20,::20],lat_grid[::20,::20], final_adjusted_U[::20,::20], final_adjusted_V[::20,::20], latlon=True, width=0.001)
# im.axes.set_aspect(1.7)
# cbar = plt.colorbar(label='wind speed [m/s]',orientation='horizontal')
# cbar.solids.set_edgecolor('face')
# plt.title("HIGH-RESOLUTION ANALYSIS VECTOR WIND FiELD | %s | %s" %(plot_timestamp,plot_tag))
# plt.savefig(fig_subdir+'wind_'+hs_plot, format='pdf')
# print('......... Corrected vector wind fields saved as: %s' %('wind_'+hs_plot))
# plt.close()

#plot DA-corrected rain field
bm_fig = plot_domain(bm)
plt.title('HIGH-RESOLUTION PRECIP ANALYSIS | %s | %s' %(plot_timestamp,plot_tag))
im = bm.imshow(elevation.T, cmap=plt.cm.gist_yarg, alpha = 0.7, origin='upper' )
im.axes.set_aspect(1.7)
im2 = bm.imshow(final_adjusted_R_plot.T, origin='upper', alpha = 0.7, cmap=cmR, norm=LogNorm(vmin=Rmin, vmax=Rmax))
im2.axes.set_aspect(1.7)
im2.set_clim(Rmin,Rmax)
cbar = plt.colorbar(im2,label='liquid precip [mm]',orientation='horizontal', ticks = rticks)
cbar.ax.set_xticklabels(rticks[1:])
scat = bm.scatter(obsTrainR['x'],obsTrainR['y'],linewidth='0.5',alpha = 0.7, s=17,marker='o', c=np.array(obsTrainR['r']), cmap=cmR, norm=LogNorm(vmin=Rmin, vmax=Rmax))
scat.set_clim(Rmin,Rmax)
plt.savefig(fig_subdir+'rain_'+hs_plot, format='pdf')
print('......... Corrected rainfall fields saved as: %s' %('rain_'+hs_plot))
print('Run COMPLETE.')

# NOTES FOR FUTURE OPTIMIZATION
# -lapse_adjust - avoid point by point (work with arrays instead)
