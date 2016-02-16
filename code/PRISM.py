# import gdal
import numpy as np
import matplotlib.pyplot as plt
from Scientific.IO import NetCDF
from scipy.spatial import KDTree
from scipy import stats
from scipy.spatial import distance
from operator import itemgetter
import os.path
from collections import Counter

# gdal.UseExceptions()
# dem = gdal.Open('cdem_aspect_151112_151410.tif')
# band = dem.GetRasterBand(1)
# aspect = band.ReadAsArray()

# #------from da_congif--------
# lat = [48.2,51.]				#[min lat, max lat]
# lon = [-128.8,-121.]			#[min lon, max lon]
# rain_roi = 0.1
# #------end--------


prism_args = np.load('prism_args.npy').item()
lon_grid, lat_grid, obs, lat, lon, rain_roi, missing_ids = itemgetter('lon_grid','lat_grid','obs','lat','lon','rain_roi','missing_ids')(prism_args)

npy_dir = './npy/'
weights_filename = 'prism_weights_%s_%s_%s_%s_r%s.npy' %(lon[0],lon[1],lat[0],lat[1],rain_roi)
lcn_filename = 'prism_lcn_%s_%s_%s_%s_r%s.npy' %(lon[0],lon[1],lat[0],lat[1],rain_roi)
id_filename = 'prism_id_%s_%s_%s_%s_r%s.npy' %(lon[0],lon[1],lat[0],lat[1],rain_roi)

#load prism data, create latlon grid arrays
nc_data = NetCDF.NetCDFFile('pr_monClim_PRISM_historical_run1_198101-201012.nc.nc', 'r')  

#subset prism data to domain bounds
min_lat = np.argmin(abs(nc_data.variables['lat'][:] - lat[0]))
max_lat = np.argmin(abs(nc_data.variables['lat'][:] - lat[1]))
min_lon = np.argmin(abs(nc_data.variables['lon'][:] - lon[0]))
max_lon = np.argmin(abs(nc_data.variables['lon'][:] - lon[1]))

dimlat = abs(max_lat - min_lat)
dimlon = abs(max_lon - min_lon)
lat_prism = np.array([nc_data.variables['lat'][min_lat:max_lat]] * dimlon)
lon_prism = np.array([nc_data.variables['lon'][min_lon:max_lon]] * dimlat).T

prism_data = nc_data.variables['pr'][:-1,min_lat:max_lat,min_lon:max_lon]
temp_mask = np.mean(prism_data,0).T
lon_prism_flat, lat_prism_flat = lon_prism[temp_mask>0], lat_prism[temp_mask>0]


#get current domain data and observatios
for key in obs:
	obs[key] = obs[key][missing_ids]

dimx, dimy = np.shape(lon_grid)
obs_lcn = zip(obs['x'],obs['y'])
prism_lcn = zip(lon_prism_flat,lat_prism_flat)
dem_lcn = zip(lon_grid.ravel(),lat_grid.ravel())

#build two KD trees for dem and prism
num_stn = len(obs_lcn)
demTree= KDTree(dem_lcn)
prismTree= KDTree(prism_lcn)
roi_dem_lcn = demTree.query_ball_point(obs_lcn, rain_roi)
roi_prism_lcn = prismTree.query_ball_point(obs_lcn, rain_roi)

#figure out the proper roi size (exclude boundary effects)`
grid_cnt = Counter([len(item) for item in roi_dem_lcn])
ballpt_size = grid_cnt.most_common(1)[0][0] 				
prism_weights = np.empty((num_stn,12,ballpt_size)) * np.nan

print('Initializing PRISM-DA routine')
for nMonth in range(12):
	print('......... Processing month %s' %(nMonth+1))
	prism_flat = prism_data[nMonth,:,:].T
	prism_flat = prism_flat[temp_mask>0]
	for nStn in range(num_stn):
		print('................. station %s out of %s' %((nStn+1), num_stn))
		if len(roi_dem_lcn[nStn]) == ballpt_size:
			select_prism_lcn = np.array(prism_lcn)[roi_prism_lcn[nStn]]									#get locations of prism and dem points surrounding obs station
			select_dem_lcn = np.array(dem_lcn)[roi_dem_lcn[nStn]]
			norm_prism_weights = prism_flat[roi_prism_lcn[nStn]]/prism_flat[roi_prism_lcn[nStn][0]]		#normalize prism weights to value at obs
			miniTree = KDTree(select_prism_lcn) 														#find nearest dem neighbour for each prism point
			prism_dist, prism_id = miniTree.query(select_dem_lcn,k=1)
			obs_coord = obs_lcn[nStn]
			for nGrid,grid in enumerate(roi_dem_lcn[nStn]): 
				pt_coord = dem_lcn[grid]
				dem_dist = distance.euclidean(pt_coord,obs_coord)
				Wdist = ((rain_roi-dem_dist)/rain_roi)**3
				Wprism = norm_prism_weights[prism_id[nGrid]]
				Wtotal = Wdist * Wprism
				prism_weights[nStn,nMonth,nGrid] = Wtotal
		else:
			print('...................... excluded due to proximity to data boundary')

#empty out location data for all incomplete roi's, force identical size
for nStn in range(num_stn):
	if len(roi_dem_lcn[nStn]) != ballpt_size:
		roi_dem_lcn[nStn] = np.empty((ballpt_size))*np.nan


# #plot newly added station affects
# plt.figure()
# plt.title('Station Locations and MD Weights', fontweight='bold')
# for nStn in range(num_stn):
# 	temp = np.reshape(prism_weights[nStn,0,:], np.shape(lon_grid))
# 	plt.imshow(temp.T, origin='upper', alpha = 0.5, cmap=plt.cm.cubehelix_r)
# cbar = plt.colorbar(label='MD weight',orientation='horizontal')
# cbar.solids.set_edgecolor('face')
# # plt.savefig('md_weights_current_run.pdf', format='pdf')
# # plt.close()
# plt.show()
		
#CURRENT SAVING SCHEME:
#if a weights file already exists for current configuration append the dictionary with new stations
weights_path = npy_dir + weights_filename
lcn_path = npy_dir + lcn_filename
id_path = npy_dir + id_filename

if os.path.isfile(weights_path) and os.path.isfile(id_path) and os.path.isfile(lcn_path):
	print '......... Appending existing prism weights file'
	prism_weights_old = np.load(weights_path)
	prism_id_old = np.load(id_path)
	prism_lcns_old = np.load(lcn_path)

	prism_weights = np.concatenate((prism_weights_old,prism_weights))
	prism_id = np.concatenate((prism_id_old,obs['id'][:]))
	prism_lcns = np.concatenate((prism_lcns_old,roi_dem_lcn))

else:
	prism_id = obs['id'][:]
	prism_lcns = roi_dem_lcn
np.save(weights_path, prism_weights)
np.save(lcn_path, prism_lcns)
np.save(id_path, prism_id)

print('......... PRISM weights saved in: ' + weights_path)
