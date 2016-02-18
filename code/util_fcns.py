def import_wrf_data(wrf_data,vars_3d,vars_4d):
	'''Import data from WRF output file (NetCDF). Assumes instantaneous values.
	(i.e. removes the first 'time' dimension)

	Parameters:
	-wrf_data: string
		path to wrf data file
	-vars_3d: tuple of strings
		string names of 4D variables
	-vars_4d: tuple of strings
		string names of 4D variables

	Returns: 
	-fc_data: dictionary 
		dictionary containing all requested variables (var names as keys)

	may 2015, nmoisseeva@eos.ubc.ca
	'''
	from Scientific.IO import NetCDF
	import numpy as np

	print('Extracting NetCDF data from: %s ' %wrf_data)
	nc_data = NetCDF.NetCDFFile(wrf_data, 'r')             
	fc_data = {}

	#copy variables of interest into a dictionary
	for var3 in vars_3d:
		print('......... %s' %var3)
		fc_data[var3] = nc_data.variables[var3][0,:,:]
	for var4 in vars_4d:
		print('......... %s' %var4)
		fc_data[var4] = nc_data.variables[var4][0,:,:,:]

	#if there are accumulated precip values included, load and subtract info from previous hour
	if ('RAINNC' in vars_3d or 'SNOWNC' in vars_3d):
		hr_back = int(wrf_data[-8:-6]) - 1
		old_wrf_data = wrf_data[:-8] + '%.2d' %hr_back + wrf_data[-6:]
		old_nc_data = NetCDF.NetCDFFile(old_wrf_data, 'r')
		fc_data['RAINNC'] = fc_data['RAINNC'] - old_nc_data.variables['RAINNC'][0,:,:]
		fc_data['SNOWNC'] = fc_data['SNOWNC'] - old_nc_data.variables['SNOWNC'][0,:,:]
		print('......... RAINNC and SNOWNC adjusted to 1H')
	
	# #if wind values included, rotate to earth coordinates:
	# if ('U10' in vars_3d or 'V10' in vars_3d):
	# 	stand_lon = float(nc_data.STAND_LON)
	# 	if 'XLONG' not in vars_3d:
	# 		print('ERROR: XLONG variable is required for wind rotation')
	# 	diff = fc_data['XLONG'] - stand_lon
	# 	if diff > 180.:
	# 		diff = diff - 360.
	# 	elif diff < -180.:
	# 		diff = diff + 360.
	# 	alpha = diff * np.pi / 180.
		
	# 	u10 = fc_data['V10']*np.sin(alpha) + fc_data['U10'] * np.cos(alpha)
	# 	v10 = fc_data['V10']*np.cos(alpha) - fc_data['U10'] * np.sin(alpha)

	# 	fc_data['U10'][:,:] = u10
	# 	fc_data['V10'][:,:] = v10
	# 	print('......... U10 and V10 rotated to earth coordinates')


	return fc_data

def import_dem_data(geotiff,lat,lon):
	'''Import dem data from geotiff file.

	Parameters:
	-geotiff: string
		path to geoTIFF data file
	-lat: list
		-[min lat, max lat]
	-lon: list
		-[min lon, max lon]

	Returns: 
	-tiffdata: array (2D)
		dem data as numpy array
	-lat_grid, lon_grid: array (2D)
		geographic coordinate arrays of elevation points 
	-bounds: tuple
		actual grid boundaries

	may 2015, nmoisseeva@eos.ubc.ca
	'''

	import gdal
	import numpy as np


	gdal.UseExceptions()
	dem = gdal.Open(geotiff)
	band = dem.GetRasterBand(1)
	tiffdata = band.ReadAsArray()

	#create lat/lon grids
	gt = dem.GetGeoTransform()							#get raster data			
	nPy, nPx = np.shape(tiffdata)						#get size of data
	lat_grid = np.zeros((nPx,nPy))						#create storage arrays
	lon_grid = np.zeros((nPx,nPy))
	for nY in range(nPy):
		for nX in range(nPx):
			lat_grid[nX,nY] = gt[3] + nX*gt[4] + nY*gt[5]
			lon_grid[nX,nY] = gt[0] + nX*gt[1] + nY*gt[2]

	min_lat_idx = np.argmin(np.abs(lat_grid[0,:]-lat[0])) 
	max_lat_idx = np.argmin(np.abs(lat_grid[0,:]-lat[1]))
	min_lon_idx = np.argmin(np.abs(lon_grid[:,0]-lon[0]))
	max_lon_idx = np.argmin(np.abs(lon_grid[:,0]-lon[1]))

	bounds = (lon_grid[min_lon_idx,0],lon_grid[max_lon_idx,0], lat_grid[0,min_lat_idx],lat_grid[0,max_lat_idx])

	if any((min_lat_idx,max_lat_idx)) == any((0,nPy)):
		print("WARNING: requestested latitude range extends to the edge (or beyond) the available geoTIFF domain")
	if any((min_lat_idx,max_lat_idx)) == any((0,nPx)):
		print("WARNING: requestested longitude range extends to the edge (or beyond) the available geoTIFF domain")

	minY, maxY = min(max_lat_idx,min_lat_idx), max(max_lat_idx,min_lat_idx)
	minX, maxX = min(max_lon_idx,min_lon_idx), max(max_lon_idx,min_lon_idx)
	lat_grid = lat_grid[minX:maxX, minY:maxY]
	lon_grid = lon_grid[minX:maxX, minY:maxY]
	tiffdata = tiffdata[minY:maxY, minX:maxX]

	tiffdata  = tiffdata.T

	return tiffdata, lat_grid, lon_grid, bounds

def prep_da(obs,lon_grid,lat_grid,elevation):
	'''Construct KDTree from dem points for future use. Remove obs stations with large elevation discrepancy 

	Parameters:
	-obs: dictionary
		obs data dictionary (from import_emwxnet_data)
	-lon_grid,lat_grid: array (2D)
		dem grid point locations (in geographic coordinates)


	Returns: 
	-obs_clean: dictionary
		cleaned obs data dictionary (from import_emwxnet_data)
	-demTree: KDTree object
		KDTree of dem points

	july 2015, nmoisseeva@eos.ubc.ca
	'''

	from scipy.spatial import KDTree
	import numpy as np

	demPoints = zip(lon_grid.ravel(),lat_grid.ravel())
	demTree = KDTree(demPoints)

	test_lcn = zip(obs['x'],obs['y'])
	dem_dist, dem_id = demTree.query(test_lcn)

	#test for bad points (elevation difference more than a set threshold)
	obs_len = len(test_lcn)
	good_pts, good_ids = [],[]
	for nStn in range(obs_len):
		idxDEM = dem_id[nStn]
		del_h = elevation.ravel()[idxDEM] - obs['h'][nStn]
		if abs(del_h) < 200:
			good_pts.append(nStn)
			good_ids.append(idxDEM)
		else:
			print('......... Excluded station ID %d: elevation mismatch = %s m [DEM: %.2f, obs: %.2f]' %(obs['id'][nStn], del_h, elevation.ravel()[idxDEM],obs['h'][nStn]))

	bad_pts = obs_len - len(good_pts)
	obs_clean = obs
	for key, val in obs.iteritems():
		obs_clean[key] = val[good_pts]
	obs_clean['dem_idx'] = good_ids

	print('......... total excluded: %s' %bad_pts)

	return obs_clean, demTree

def verif_sets(fcPoints,demTree,obs,fcVAR,demVAR,tag,plot_timestamp,*args,**kwarg):

	'''Perform point-by-point verification of forecast and observations. Compare
	to downscaled forecast and observations.

	Parameters:
	-fcPoints: tuple (lon, lat)
		forecast location points (in geographic coordinates)
	-demTree: KDTree object
		KDtree object based on lon/lat grid points
	-obs: dictionary
		obs data dictionary (from import_emwxnet_data)
	-fcVAR: array (2D)
		forecast variable 
	-demVAR: array (2D)
		downscaled forecast variable 
	-tag: string
		method being tested (for reference/title only)
	-plot_timestamp: string
		wrf output date and time 
	-land_test: boolean
		flag to use land stations only for verification; requires dem_landmask argument (2D, high-resolution)
	-var: string
			dictionary key for corresponding variable

	Returns: 
	-verif_fig: figure
		scatter subplots of fc vs obs and downscaled fc vs obs 

	may 2015, nmoisseeva@eos.ubc.ca
	'''
	from scipy.spatial import KDTree
	from scipy import stats
	import numpy as np
	from matplotlib import cm
	import matplotlib.pyplot as plt

	key = kwarg['var']
	test_lcn = zip(obs['x'],obs['y'])
	dem_dist, dem_id = demTree.query(test_lcn)

	#if set in LAND_TEST mode exclude all on-water points based on dem_landmask
	if kwarg['land_test']:
		print('WARNING: verification set in LAND_TEST mode - all on-water locations will be excluded')
		dem_landmask = args[0]
		landmask_flat = dem_landmask.ravel()
		land_idx = np.where(landmask_flat[dem_id]==1)[0]
		test_lcn = np.array(test_lcn)[land_idx]
		print('......... total number of land points used for verification: %s' %len(test_lcn))
		dem_dist, dem_id = demTree.query(test_lcn)

	#construct second KD-tree for forecast points
	fcTree = KDTree(fcPoints)
	fc_dist,fc_id = fcTree.query(test_lcn)
	obs_len = len(test_lcn)

	ptFCval, ptDEMval = [],[]
	fcVAR_flat, demVAR_flat  = fcVAR.ravel(), demVAR.ravel()

	#get values closest to obs point from forecast and high-rest dataset
	for nPt in range(obs_len):
		idxFC, idxDEM = fc_id[nPt], dem_id[nPt]
		ptFCval.append(fcVAR_flat[idxFC])
		ptDEMval.append(demVAR_flat[idxDEM])

	ptOBSval = np.asarray(obs[key])			#get values at obs points
	h_cmap = obs['h']						#save elevation for colormap
	if kwarg['land_test']:
		ptOBSval = ptOBSval[land_idx]
		h_cmap = np.array(h_cmap)[land_idx]
	ptFCval, ptDEMval = np.asarray(ptFCval), np.asarray(ptDEMval)

	#perform some stats
	rmseFC = np.sqrt(np.nanmean((ptOBSval - ptFCval)**2))
	rmseDEM = np.sqrt(np.nanmean((ptOBSval - ptDEMval)**2))
	mf,bf,rf,pf,stdf = stats.linregress(ptOBSval,ptFCval)
	md,bd,rd,pd,stdd = stats.linregress(ptOBSval,ptDEMval)
	maeFC = np.nanmean(abs(ptOBSval - ptFCval))
	maeDEM = np.nanmean(abs(ptOBSval - ptDEMval))
	biasFC = np.nanmean(ptOBSval) - np.nanmean(ptFCval)
	biasDEM = np.nanmean(ptOBSval) - np.nanmean(ptDEMval)

	#plotting
	verif_fig = plt.figure(figsize=(20, 10))
	plt.subplot(1,2,1)
	plt.suptitle('\n Verification of %s: %s' %(tag,plot_timestamp), fontsize=10)
	plt.title('Observations vs Raw Forecast', fontsize=10, fontweight='bold')
	lims = (min(min(ptOBSval),min(ptDEMval),min(ptFCval)),max(max(ptOBSval),max(ptDEMval),max(ptFCval)))
	plt.scatter(ptFCval, ptOBSval, linewidth='0', c=h_cmap, cmap=plt.cm.coolwarm)
	plt.xlabel('raw forecast')
	plt.ylabel('observations')
	plt.xlim(lims)
	plt.ylim(lims)
	plt.annotate('\n RMSE: %.2f \n\n R: %.2f \n\n MAE: %.2f \n\n Bias: %.2f \n' %(rmseFC,rf,maeFC,biasFC), xy=(0.05, 0.72), xycoords='axes fraction', \
		bbox=dict(facecolor='none', edgecolor='black'), fontsize=14, fontweight='bold')
	plt.subplot(1,2,2)
	plt.title('Observations vs Adjusted Forecast', fontsize=10, fontweight='bold')
	scat = plt.scatter(ptDEMval, ptOBSval, linewidth='0', c=h_cmap, cmap=plt.cm.coolwarm)
	plt.xlabel('adjusted forecast')
	plt.ylabel('observations')
	plt.xlim(lims)
	plt.ylim(lims)
	plt.annotate('\n RMSE: %.2f  \n\n R: %.2f \n\n MAE: %.2f \n\n Bias: %.2f \n' %(rmseDEM,rd,maeDEM,biasDEM), xy=(0.05, 0.72), xycoords='axes fraction', \
		bbox=dict(facecolor='none', edgecolor='black'), fontsize=14, fontweight='bold')
	fig = plt.gcf()
	cax = fig.add_axes([0.93, 0.1, 0.01, 0.8])
	cbar = fig.colorbar(scat, cax=cax)
	cbar.set_label('elevation [m]', rotation=90)

	return verif_fig

def import_emwxnet_data(emx_name,emx_dir,bounds,hr):
	'''Import observation data from EmWxNet.

	Parameters:
	-emx_name: string
		name of file containing list of stations and locations
	-emx_dir: string
		path to EmWxNet data folder
	-bounds: tuple
		domain boundaries (min lon, max lon, min lat, max lat)
	-hr: string
		hour for which to extract observations 

	Returns: 
	-obsT: dictionary
		{'x': longitudes, 'y': latitudes, 't': temperatures, 'h': heights, 'id': station ids}
	-obsT_len: int
		number of stations with usable temperature data
	-obsW: dictionary
		{'x': longitudes, 'y': latitudes, 'h': heights, 'id': station ids,'u': wind component, 'v': wind component}
	-obsW_len: int
		number of stations with usable wind data

	october 2015, nmoisseeva@eos.ubc.ca
	'''

	import numpy as np
	import os.path
	from datetime import datetime
	import matplotlib.path as mplPath

	obsfile = os.path.join(emx_dir,emx_name)

	print('Extracting observation locations from: %s' %obsfile)
	obs_lcn = np.genfromtxt(obsfile, usecols=(0,1,2,3), skip_header=1, autostrip=True,dtype=float)

	num_stn = np.shape(obs_lcn)[0]						#get number of stations
	time0 = datetime.strptime(hr,'%H')

	print('......... total number of stations: %s' %num_stn)
	# obs_id, obs_x, obs_y, obs_h, obs_t, obs_r = [],[],[],[],[],[]						#storage arrays for location coordinates

	obsT = {'x':np.array([]),'y': np.array([]), 't': np.array([]), 'h':np.array([]), 'id':np.array([])}
	obsW = {'x':np.array([]),'y': np.array([]), 'u': np.array([]), 'v': np.array([]), 'h':np.array([]), 'id':np.array([])}
	obsR = {'x':np.array([]),'y': np.array([]), 'r': np.array([]), 'h':np.array([]), 'id':np.array([])}


	#loop through all available observation stations, store stns with data and in Canada only
	for nStn in range(num_stn):	
		stn_file = os.path.join(emx_dir,str(int(obs_lcn[nStn,0]))+'.txt')
		stnID, stnLat, stnLon, stnH = obs_lcn[nStn,:]
		#make sure station exists and is not in the US territory, or outside of domain
		if os.path.isfile(stn_file) \
			and not(stnLon>-123.03 and stnLat<49.) \
			and (stnLon>bounds[0] and stnLon<bounds[1]) \
			and (stnLat>bounds[2] and stnLat<bounds[3]):
			f = open(stn_file)
			row_count = len(f.readlines())
			if row_count > 2:													#make sure file is not empty
				obs_data = np.genfromtxt(stn_file, \
					dtype = [('date',int),('time','S6'),('temp',float),('wind_speed',float),('wind_dir',float),('rain',float)], skip_header=1)
				obs_time_tag = [x.zfill(6) for x in obs_data['time'][:]]		#add leading zeros to format as HHMMSS
				timeX = [datetime.strptime(x, '%H%M%S') for x in obs_time_tag]	#convert to datetime
				del_time = [abs((time0 - x).total_seconds()/60.) for x in timeX ]	#calculate time offset in minutes
				closest_t_idx = np.argmin(del_time)								#get index of closest time stamp
				pt_temp = obs_data['temp'][closest_t_idx]						#get temperature reading
				pt_rain = obs_data['rain'][closest_t_idx]						#get rain reading

				if (del_time[closest_t_idx] < 5):								#ensure reasonable del_t range
					if (-50 < pt_temp < 60):									#ensure reasonable temp
						obsT['x'] = np.append(obsT['x'],stnLon)
						obsT['y'] = np.append(obsT['y'],stnLat)
						obsT['h'] = np.append(obsT['h'],stnH)
						obsT['id'] = np.append(obsT['id'],stnID)
						obsT['t'] = np.append(obsT['t'],pt_temp+273.)
					if (0. <= pt_rain < 100):									#ensure reasonable precip
						obsR['x'] = np.append(obsR['x'],stnLon)
						obsR['y'] = np.append(obsR['y'],stnLat)
						obsR['h'] = np.append(obsR['h'],stnH)
						obsR['id'] = np.append(obsR['id'],stnID)
						obsR['r'] = np.append(obsR['r'],pt_rain)

				#get wind data as an average of all values within the past hour
				temp_u, temp_v, cnt = 0, 0, 0
				for nTime, timerec in enumerate(timeX):
					if timerec.hour == time0.hour:
						pt_wdir = obs_data['wind_dir'][nTime]     				#get wind direction
						pt_wspd = obs_data['wind_speed'][nTime]					#get wind speed in (ASSUMES m/s)
						if ( 6. < pt_wspd < 70.): 	 							#ensure reasonable sensitivity range
							u = (pt_wspd/3.6) * np.cos(np.radians(270-pt_wdir)) 	#convert to components
							v = (pt_wspd/3.6) * np.sin(np.radians(270-pt_wdir))
							temp_u = temp_u + u
							temp_v = temp_v + v
							cnt = cnt + 1.
							# print pt_wspd,pt_wdir, u, v
				if cnt > 0:
					obsW['x'] = np.append(obsW['x'],stnLon)
					obsW['y'] = np.append(obsW['y'],stnLat)
					obsW['h'] = np.append(obsW['h'],stnH)
					obsW['id'] = np.append(obsW['id'],stnID)
					pt_u = temp_u/cnt
					pt_v = temp_v/cnt
					obsW['u'] = np.append(obsW['u'],pt_u)
					obsW['v'] = np.append(obsW['v'],pt_v)

	obsT_len = len(obsT['x'])
	obsW_len = len(obsW['x'])
	obsR_len = len(obsR['x'])

	print('......... total stations with useful temperature data: %s' %obsT_len)
	print('......... total stations with useful wind data: %s' %obsW_len)
	print('......... total stations with useful rain data: %s' %obsR_len)
	
	return obsT, obsT_len, obsW, obsW_len, obsR, obsR_len

def plot_domain(bm):
	'''Plots entire analysis domain. Masks US territory. 

	Parameters:
	-bm: mpl basemap
		current basemap object

	Returns: 
	-bm_fig: figure
		coastal outlines and masked US 

	may 2015, nmoisseeva@eos.ubc.ca
	'''
	from mpl_toolkits import basemap
	from matplotlib.patches import Polygon
	import matplotlib.pyplot as plt


	bm_fig = plt.figure(figsize=(15, 13))
	ax = plt.gca()
	#mask US territory
	mask = us_mask(bm)
	poly = Polygon(mask, facecolor='grey',linewidth=None, edgecolor='grey')
	ax.add_patch(poly)	
	#add watersheds
	bm.readshapefile('./shape_files/MajorHydroWatershedsProject',name='watersheds',drawbounds=True, linewidth=0.4, color='#583100')
	bm.drawcoastlines(linewidth=0.3,color='grey')

	return bm_fig

def us_mask(bm):
	'''Gets verticies of US territory boundaries(can be used as mask). Saves for future use in the run. 

	Parameters:
	-bm: mpl basemap
		current basemap object

	Returns: 
	mask: list (2D)
		list of verticies of us boundaries

	june 2015, nmoisseeva@eos.ubc.ca
	'''
	import numpy as np
	from mpl_toolkits import basemap

	bm.readshapefile('./shape_files/st99_d00', name='states', drawbounds=False)
	state_names = []
	for shape_dict in bm.states_info:
	    state_names.append(shape_dict['NAME'])
	mask = bm.states[state_names.index('Washington')]

	return mask

def get_gamma(wrfT,wrfP,wrfPB,wrfPH,wrfPHB,lvl):
	'''Calculates a gamma field (2D) for each surface model grid point and 
	returns a flattend array. 

	Parameters:
	-wrfT, wrfP,wrfPB,wrfPH,pwrfPHB: array (3D)
		wrf variable arrays [height level, y-grid, x-grid]
	-lvl: int
		maximum model level to use for lapse rate calculation 

	Returns: 
	-fcGamma_flat: array (1D)
		flattened 

	may 2015, nmoisseeva@eos.ubc.ca
	'''
	import numpy as np
	from scipy import stats

	g = 9.81

	lY, lX = np.shape(wrfT)[1:3]
	fcLevelHGT = (wrfPH + wrfPHB)/g
	fc_lapse_field = np.empty((lY,lX))

	for nY in range(lY):
		for nX in range(lX):
			#convert potential temperature to "regular" temperature
			testT = (wrfT[0:lvl,nY,nX] + 300) * ((wrfP[0:lvl,nY,nX]+wrfPB[0:lvl,nY,nX])/100000)**(2./7)
			testH = fcLevelHGT[0:lvl,nY,nX] 
			m,b,r,p,std_err = stats.linregress(testT,testH)
			fc_lapse_field[nY,nX] = 1000/m
	fcGamma_flat = fc_lapse_field.ravel()

	return fcGamma_flat

def lapse_adjust(interpVAR,elevation,interpH,interpGamma):
	'''Performs lapse rate adjustment of the supplied variable field. 

	Parameters:
	-elevation: array (2D)
		high resolution DEM array
	-interpVAR, interpH, interpGamma: array (2D)
		interpolated (basic) height, lapse rate and variable of interest

	Returns: 
	-demVAR: array (2D)
		downscaled lapse-rate adjusted variable field

	may 2015, nmoisseeva@eos.ubc.ca
	'''
	import numpy as np

	interpGamma_flat = interpGamma.ravel()
	demH_flat = elevation.ravel()
	data_len = len(demH_flat)

	interpVAR_flat, interpH_flat = interpVAR.ravel(),interpH.ravel()

	demVAR_flat = np.zeros(data_len)*np.nan
	for nPt in range(data_len):
		delH = demH_flat[nPt] - interpH_flat[nPt]
		del_VAR = delH*interpGamma_flat[nPt]/1000.

		demVAR_flat[nPt] = interpVAR_flat[nPt] + del_VAR
	demVAR = np.reshape(demVAR_flat, np.shape(elevation))

	return demVAR

def subset_obs(obs, verif_frac):
	'''Subsets observations for cross-evaluation.

	Parameters:
	-obs: dictionary
		obs data dictionary (cleaned)
	-verif_frac: float
		fraction of data to use for training (set by user in da_config.py)

	Returns: 
	-obsTest, obsTrain: dictionary
		training and test data set dictionaries

	may 2015, nmoisseeva@eos.ubc.ca
	'''
	import numpy as np
	import random 

	obs_len = len(obs['x'])
	rand_idx = random.sample(np.arange(0,obs_len),int(obs_len*verif_frac))
	test_idx = list(set(np.arange(0,obs_len)) - set(rand_idx))

	obsTrain = {}
	obsTest = {}
	for key, val in obs.iteritems():
		obsTrain[key] = np.array(obs[key])[rand_idx]
		obsTest[key] = np.array(obs[key])[test_idx]

	return obsTrain, obsTest

def da_roi(roi,elev_roi,obsTrain,elevation,demVAR,lon_grid,lat_grid,dem_landmask,bias_mode,*args,**kwarg):
	'''Basic data assimilation using horizontal and vertical ROI regions. 

	Parameters:
	-roi: float
		horizontal region of influence in decimal degrees
	-elev_roi: float
		vertical region of influence in meters
	-verif_frac: float
		percentage of data to use for assimilation (the rest for verification)
	-obs: dictionary
		obs data dictionary (from import_emwxnet_data)
	-elevation: 2D array 
		-dem data from geoTIFF
	-demVAR: 2D array
		-downscaled variable
	-lon_grid, lat_grid: 2D array
		dem lon/lat grids
	-interpGamma: 2D array 
		downscaled lapse rate (for temperature assimilation only)
	-kwarg 'var': string
		diction key of corresponding variable

	Returns: 
	 -final_adjusted_T: array (2D)
	 	final high resolution field with assimilated data
	june 2015, nmoisseeva@eos.ubc.ca
	'''

	import numpy as np
	from scipy.spatial import KDTree
	from scipy import stats
	from scipy.spatial import distance

	key = kwarg['var']

	if key=='t' and bias_mode==False:
		if len(args)==0:
			print('ERROR: interpGamma missing for temperature assimilation')
			sys.exit()
		interpGamma = args[0]
		demGamma_flat = interpGamma.ravel()

	#import high resolution dem 
	demH_flat = elevation.ravel()
	demVAR_flat = demVAR.ravel()
	dem_landmask_flat = dem_landmask.ravel()

	#first pass KD-tree 
	dem_lcn = np.array(zip(lon_grid.ravel(),lat_grid.ravel()))
	obs_lcn = np.array(zip(obsTrain['x'],obsTrain['y']))
	demTree = KDTree(dem_lcn)

	#get indecies of affected areas
	roiIDs = demTree.query_ball_point(obs_lcn, roi)
	roi_idx_raw= []

	#record indecies of affected cells (excluding water for temperature)
	for nPt in range(len(obsTrain['x'])):
		idxVals = roiIDs[nPt][:]
		if key=='t':
			idxVals = np.array(idxVals)[dem_landmask_flat[idxVals].astype(bool)]
		roi_idx_raw.extend(idxVals)

	#get indecies exluding repeats from multiple sources
	roi_idx = set(roi_idx_raw)
	roi_idx = list(roi_idx)
	affected_dem = dem_lcn[roi_idx]
	affected_H = demH_flat[roi_idx]
	affected_VAR = demVAR_flat[roi_idx]
	if key=='t' and bias_mode==False:
		affected_Gamma = demGamma_flat[roi_idx]

	#second pass KDtree
	obsTree = KDTree(obs_lcn)
	neighbourIDs = obsTree.query_ball_point(affected_dem, roi)

	#perform correction
	weightedVAR = []
	Wf = 1
	if bias_mode:
		nearest_VAR = demVAR_flat[obsTrain['dem_idx'][:]]
		for nPt in range(len(neighbourIDs)):
			pairs = neighbourIDs[nPt]
			pt_coord = affected_dem[nPt]
			pt_height = affected_H[nPt]
			pt_VAR = affected_VAR[nPt]
			for nPair in range(len(pairs)):
				obs_idx = pairs[nPair]
				obs_coord = obs_lcn[obs_idx]
				obs_height = obsTrain['h'][obs_idx]
				obs_VAR = obsTrain[key][obs_idx]
				dist = distance.euclidean(pt_coord,obs_coord)
				bias = obs_VAR - nearest_VAR[obs_idx]
				Wd = ((roi-dist)/roi)**2
				vert_dist = pt_height - obs_height
				We = ((elev_roi-abs(vert_dist))/elev_roi)**2
				WT = Wd*We
				corrected_val = pt_VAR + bias * WT
				pt_VAR = corrected_val
			weightedVAR.append(pt_VAR)
	else:
		for nPt in range(len(neighbourIDs)):
			pairs = neighbourIDs[nPt]
			pt_coord = affected_dem[nPt]
			pt_height = affected_H[nPt]
			pt_VAR = affected_VAR[nPt]
			if key=='t':
				pt_lapse = affected_Gamma[nPt]
			for nPair in range(len(pairs)):
				obs_idx = pairs[nPair]
				obs_coord = obs_lcn[obs_idx]
				obs_height = obsTrain['h'][obs_idx]
				obs_VAR = obsTrain[key][obs_idx]
				dist = distance.euclidean(pt_coord,obs_coord)
				Wd = ((roi-dist)/roi)**2
				vert_dist = pt_height - obs_height
				We = ((elev_roi-abs(vert_dist))/elev_roi)**2
				WT = Wd*We
				if key=='t':
					corrected_val = pt_VAR*(1-WT) + (obs_VAR+ (vert_dist*pt_lapse/1000))*WT
				else:
					corrected_val = pt_VAR*(1-WT) + obs_VAR*WT 
				pt_VAR = corrected_val
			weightedVAR.append(pt_VAR)

	weightedVAR = np.array(weightedVAR)
	full_adjusted_VAR = np.copy(demVAR_flat)
	full_adjusted_VAR[roi_idx] = weightedVAR
	final_adjusted_VAR = np.reshape(full_adjusted_VAR, np.shape(elevation))

	return final_adjusted_VAR

def make_landmask(bm,lon_grid,lat_grid):
	'''Performs lapse rate adjustment of the supplied variable field. 

	Parameters:
	-bm: basemaps object
		basemap configuration to be used for creating a landmask
	-interpVAR, interpH, interpGamma: array (2D)
		interpolated (basic) height, lapse rate and variable of interest

	Returns: 
	-demVAR: array (2D)
		downscaled lapse-rate adjusted variable field

	may 2015, nmoisseeva@eos.ubc.ca
	'''

	from matplotlib import path 
	import numpy as np

	dem_lcn = np.array(zip(lon_grid.ravel(),lat_grid.ravel()))
	polygons = bm.coastpolygons
	poly_type = bm.coastpolygontypes

	landmask = np.zeros(len(dem_lcn))
	for nPoly in range(len(poly_type)):
		poly = polygons[nPoly]
		poly = zip(*poly)
		poly_path = path.Path(poly)
		mask = poly_path.contains_points(dem_lcn)
		if poly_type[nPoly]==1:
			landmask[mask==1] = 1
		elif poly_type[nPoly]==2:
			landmask[mask==1] = 0
	dem_landmask = np.reshape(landmask, np.shape(lon_grid))

	return dem_landmask

def split_interp(fc_data,lon_grid,lat_grid, dem_landmask,**kwarg):
	'''Performs lapse rate adjustment of the supplied variable field. 

	Parameters:
	-fc_data: dictionary
		dictionary generated from WRF output
	-lat_grid, lon_grid: array (2D)
		high resolution dem grid 
	-dem_landmask: array (2D)
		high resolution land mask 
	-var: string
		keyword corresponding to variable key (e.g. 'T2','U10') in fc_data

	Returns: 
	-interpVAR: array (2D)
		split interpolated variable field (high resolution)

	june 2015, nmoisseeva@eos.ubc.ca
	'''
	import numpy as np 
	import numpy.ma as ma
	from scipy.interpolate import griddata

	# global fcx, fcy, fcPoints, fcT, fcT_flat
	key = kwarg['var']

	fcx ,fcy = (fc_data['XLONG'][:,:], fc_data['XLAT'][:,:])
	fcPoints = zip(fcx.ravel(), fcy.ravel())
	fcH, fcVAR = fc_data['HGT'][:,:], fc_data[key][:,:]
	fcVAR_flat, fcH_flat  = fcVAR.ravel(), fcH.ravel()
	fcLand = fc_data['LANDMASK'][:,:].ravel()

	# separate water and land points to avoid edge artifacts
	fc_land_idx = np.where(fcLand==1)[0]
	fc_water_idx = np.where(fcLand==0)[0]
	fcPoints_land = np.array(fcPoints)[fc_land_idx]
	fcPoints_water = np.array(fcPoints)[fc_water_idx]
	fcVAR_land = np.array(fcVAR_flat)[fc_land_idx]
	fcVAR_water = np.array(fcVAR_flat)[fc_water_idx]

	#basic cubic interpolation of topography and temperature (no lapse rate effects)
	interpVAR_land = griddata(fcPoints_land, fcVAR_land, (lon_grid, lat_grid), method='cubic')
	interpVAR_water = griddata(fcPoints_water, fcVAR_water, (lon_grid, lat_grid), method='cubic')

	ma_water = ma.masked_array(interpVAR_water, dem_landmask)
	ma_land = ma.masked_array(interpVAR_land, np.logical_not(dem_landmask))
	interpVAR = ma_water.filled(fill_value=0) + ma_land.filled(fill_value=0)

	#if working with wind data - apply smoother on the land boundaries
	if (key=='U10' or key=='V10'):
		from scipy.signal import convolve2d
		kernel = np.ones((9,9))/81
		temp = convolve2d(interpVAR, kernel)
		interpVAR = temp[4:-4,4:-4]

	return interpVAR

def da_md(obsTrain,elevation,lon_grid,lat_grid,dem_landmask,demVAR,params,bias_mode,dist_cutoff,bm,*args,**kwarg):
	'''Mother-Daugher appreach for assimilating observations. 

	Parameters:
	-obsTrain: dictionary
		obs data dictionary
	-elevation: 2D array 
		-dem data from geoTIFF 
	-lon_grid, lat_grid: 2D array
		dem lon/lat grids
	-dem_landmask: array (2D)
		high resolution land mask 
	-demVAR: 2D array
		-downscaled variable
	-params: list
		[a,b, Z_ref1, Z_ref2] - parameters for sharing factor calculation
	-bias_mode: boolean
		1: run in bias mode (correct increment), 0: correct actual variable
	-dist_cutoff: float
		maximum anisotropic horizontal distance in degrees to continue iteration
	-bm: basemaps object
		current basemap
	-interpGamma: 2D array 
		downscaled lapse rate (for temperature assimilation only)
	-kwarg 'var': string
		diction key of corresponding variable
	-kwarg 'land_mode': boolean
		use land mode (temperature) or all grids (wind)

	Returns: 
	 -final_adjusted_VAR: array (2D)
	 	final high resolution field with assimilated data


	october 2015, nmoisseeva@eos.ubc.ca
	'''

	import numpy as np
	import pickle
	import matplotlib.pyplot as plt
	import os.path
	import sys

	file_tag = params + (dist_cutoff,)

	key = kwarg['var']
	land_mode = kwarg['land_mode']

	if key=='t':
		if len(args)==0:
			print('ERROR: interpGamma missing for temperature assimilation')
			sys.exit()

	weights_path = './npy/md_weights_%s-%s-%s-%s_r%s.npy' %file_tag
	dist_path = './npy/md_dist_%s-%s-%s-%s_r%s.npy' %file_tag
	id_path = './npy/md_id_%s-%s-%s-%s_r%s.npy' %file_tag
	obs_len = len(obsTrain['x'])
	landmask_flat = dem_landmask.ravel()
	VAR_flat = demVAR.ravel()

	#load or calculate weights
	if os.path.isfile(weights_path)==False or os.path.isfile(id_path)==False or os.path.isfile(dist_path)==False:
		print('WARNING: no existing weights/IDs found for given parameters: calculating new weights')
		save_data = {'elevation':elevation,'lon_grid':lon_grid,'lat_grid':lat_grid,'obs':obsTrain,'landmask_flat':landmask_flat,'missing_ids':np.arange(obs_len)}
		np.save('md_args.npy',save_data)	
		os.system('python2.7 MD.py')									#calculate new weights for all stations

	else:
		print('MD weights, distances and IDs found at: %s, %s, %s' %(weights_path,dist_path,id_path))
		md_id = np.load(id_path)
		if set(obsTrain['id'][:]).issubset(set(md_id))==False:
			missing_ids = [i for i,item in enumerate(obsTrain['id'][:]) if item not in md_id]		#find stations that are missing weights
			print('WARNING: Weights missing for %s station ids - updating weights file' %len(missing_ids))
			save_data = {'elevation':elevation,'lon_grid':lon_grid,'lat_grid':lat_grid,'obs':obsTrain,'landmask_flat':landmask_flat,'missing_ids':missing_ids}
			np.save('md_args.npy',save_data)
			os.system('python2.7 MD.py')								#run MD for missing stations only

	if os.path.isfile('md_args.npy'): 								#if MD was ran, clean up 
		os.system('rm md_args.npy')

	#FOR IMPLEMENTATION WITH PYTHON 3.4
	# md_dict = np.load(weights_path).item()
	# md_weights, md_id = md_dict['weights'], md_dict['id']

	md_weights = np.load(weights_path) 								#load fresh weights	
	md_dist = np.load(dist_path)									#load fresh distances
	md_id = np.load(id_path) 										#load fresh ID list

	stn_idx = [np.where(i==md_id)[0][0] for i in obsTrain['id']] 	#get ID's of stations for assimilation
	select_weights = md_weights[stn_idx,:] 							#get weights for the selected stations
	select_dist = md_dist[stn_idx,:]

	md_fix = (((dist_cutoff - select_dist)/dist_cutoff)**2) * select_weights

	#plot newly added station affects
	bm_fig = plot_domain(bm)
	bm.imshow(elevation.T, cmap=plt.cm.gist_yarg, alpha = 0.5, origin='upper' )
	bm.drawcoastlines()
	plt.title('Station Locations and MD Weights', fontweight='bold')
	for nStn in range(len(stn_idx)):
		temp = np.reshape(md_fix[nStn,:], np.shape(elevation))
		bm.imshow(temp.T, origin='upper', alpha = 0.5, cmap=plt.cm.cubehelix_r)
	cbar = plt.colorbar(label='MD weight',orientation='horizontal')
	cbar.solids.set_edgecolor('face')
	bm.scatter(obsTrain['x'],obsTrain['y'],s=6,marker='o', color ='k')
	plt.savefig('md_weights_current_run.pdf', format='pdf')
	plt.close()

	#calculate correction  and plot current run
	md_VAR = np.empty_like(select_weights) * np.nan
	Wf = 1. 											# assigned forecast weight fraction
	if bias_mode:
		for nStn in range(obs_len):
			correction_factor = md_fix[nStn,:]
			if land_mode:
				correction_factor[landmask_flat==0] = np.nan
			bias = obsTrain[key][nStn] - VAR_flat[obsTrain['dem_idx'][nStn]]
			md_VAR[nStn,:] = (np.array(VAR_flat) + bias) * correction_factor
	else:
		for nStn in range(obs_len):
			correction_factor = md_fix[nStn,:]
			if land_mode:
				correction_factor[landmask_flat==0] = np.nan
			if key=='t':
				interpGamma = args[0]
				delH = elevation.ravel() - obsTrain['h'][nStn]
				delVAR = delH *  interpGamma.ravel()/1000.
				rawVAR = obsTrain[key][nStn] + delVAR
			else:
				rawVAR = np.ones(np.shape(VAR_flat)) * obsTrain[key][nStn]
			md_VAR[nStn,:] = rawVAR[:] * correction_factor

	#combine corrections from all stations: VAR = (Wf*VARf + SUM(Wmd*VARobs)) / (Wf + SUM(Wmd))
	VARtotal = ( Wf*VAR_flat + np.nansum(md_VAR,axis=0) )/ ( Wf + np.nansum(md_fix,axis=0) )
	final_adjusted_VAR = np.reshape(VARtotal, np.shape(lon_grid))

	return final_adjusted_VAR

def da_prism(obsTrain,elevation,lon_grid,lat_grid,demVAR,month,lat,lon,rain_roi,bm,bias_mode):
	'''PRISM-corrected roi approach for assimilating observations. 

	Parameters:
	-obsTrain: dictionary
		obs data dictionary
	-elevation: 2D array 
		-dem data from geoTIFF 
	-lon_grid, lat_grid: 2D array
		dem lon/lat grids
	-demVAR: 2D array
		-downscaled variable
	-month: int
		current month
	-lat, lon: 2-element lists
		domain bounds from da_config
	-rain_roi: float
		region of influence extent from da_config
	-bias_mode: boolean
		1: run in bias mode (correct increment), 0: correct actual variable
	-bm: basemaps object
		current basemap

	Returns: 
	 -final_adjusted_VAR: array (2D)
	 	final high resolution field with assimilated data


	november 2015, nmoisseeva@eos.ubc.ca
	'''

	import numpy as np
	import pickle
	import matplotlib.pyplot as plt
	import os.path
	import sys

	weights_path = './npy/prism_weights_%s_%s_%s_%s_r%s.npy' %(lon[0],lon[1],lat[0],lat[1],rain_roi)
	lcn_path = './npy/prism_lcn_%s_%s_%s_%s_r%s.npy' %(lon[0],lon[1],lat[0],lat[1],rain_roi)
	id_path = './npy/prism_id_%s_%s_%s_%s_r%s.npy' %(lon[0],lon[1],lat[0],lat[1],rain_roi)
	obs_len = len(obsTrain['x'])
	VAR_flat = demVAR.ravel()
	dem_lcn = zip(lon_grid.ravel(),lat_grid.ravel())

	#load or calculate weights
	if os.path.isfile(weights_path)==False or os.path.isfile(id_path)==False or os.path.isfile(lcn_path)==False:
		print('WARNING: no existing PRISM weights/IDs found for given parameters: calculating new weights')
		save_data = {'lon_grid': lon_grid,'lat_grid':lat_grid,'obs':obsTrain,'lat':lat,'lon':lon,'rain_roi':rain_roi,'missing_ids':np.arange(obs_len)}
		np.save('prism_args.npy',save_data)	
		os.system('python2.7 PRISM.py')									#calculate new weights for all stations

	else:
		print('PRISM weights, locations and IDs found at: %s, %s, %s' %(weights_path,lcn_path,id_path))
		prism_id = np.load(id_path)
		if set(obsTrain['id'][:]).issubset(set(prism_id))==False:
			missing_ids = [i for i,item in enumerate(obsTrain['id'][:]) if item not in prism_id]		#find stations that are missing weights
			print('WARNING: Weights missing for %s station ids - updating weights file' %len(missing_ids))
			save_data = {'lon_grid': lon_grid,'lat_grid':lat_grid,'obs':obsTrain,'lat':lat,'lon':lon,'rain_roi':rain_roi,'missing_ids':missing_ids}
			np.save('prism_args.npy',save_data)
			os.system('python2.7 PRISM.py')								#run PRISM for missing stations only

	if os.path.isfile('prism_args.npy'): 								#if PRISM was ran, clean up 
		os.system('rm prism_args.npy')

	#FOR IMPLEMENTATION WITH PYTHON 3.4
	# md_dict = np.load(weights_path).item()
	# md_weights, md_id = md_dict['weights'], md_dict['id']

	prism_weights = np.load(weights_path) 								#load fresh weights	
	prism_lcn = np.load(lcn_path)
	prism_id = np.load(id_path) 										#load fresh ID list

	stn_idx = [np.where(i==prism_id)[0][0] for i in obsTrain['id']] 	#get ID's of stations for assimilation
	select_weights = prism_weights[stn_idx,month-1,:] 					#get weights for the selected stations
	select_lcns = prism_lcn[stn_idx]

	#calculate correction  and plot current run
	prism_VAR = np.empty((obs_len,len(VAR_flat)))* np.nan
	array_PRISM = np.empty((obs_len,len(VAR_flat)))* np.nan

	Wf = 0.5  											# assigned forecast weight fraction
	if bias_mode:
		for nStn in range(obs_len):
			correction_factor = select_weights[nStn,:]
			bias = obsTrain['r'][nStn] - VAR_flat[obsTrain['dem_idx'][nStn]]
			if np.isfinite(sum(select_lcns[nStn])):
				prism_VAR[nStn,select_lcns[nStn]] = (np.array(VAR_flat[select_lcns[nStn]]) + bias) * correction_factor
				array_PRISM[nStn,select_lcns[nStn]] = select_weights[nStn]
	else:
		for nStn in range(obs_len):
			correction_factor = select_weights[nStn,:]
			rawVAR = np.ones(np.shape(correction_factor)) * obsTrain['r'][nStn]
			if np.isfinite(sum(select_lcns[nStn])):
				array_PRISM[nStn,select_lcns[nStn]] = select_weights[nStn]
				prism_VAR[nStn,select_lcns[nStn]] = rawVAR[:] * correction_factor

	#combine corrections from all stations: VAR = (Wf*VARf + SUM(Wmd*VARobs)) / (Wf + SUM(Wmd))
	VARtotal = ( Wf*VAR_flat + np.nansum(prism_VAR,axis=0) )/ ( Wf + np.nansum(array_PRISM,axis=0) )
	final_adjusted_VAR = np.reshape(VARtotal, np.shape(lon_grid))

	#plot station affects
	bm_fig = plot_domain(bm)
	bm.imshow(elevation.T, cmap=plt.cm.gist_yarg, alpha = 0.5, origin='upper' )
	bm.drawcoastlines()
	plt.title('Station Locations and PRISM Weights', fontweight='bold')
	for nStn in range(len(stn_idx)):
		temp = np.reshape(array_PRISM[nStn,:], np.shape(elevation))
		bm.imshow(temp.T, origin='upper', alpha = 0.5)
	cbar = plt.colorbar(label='PRISM weight',orientation='horizontal')
	bm.scatter(obsTrain['x'],obsTrain['y'],s=6,marker='o', color ='k')
	plt.savefig('prism_weights_current_run.pdf', format='pdf')
	plt.close()


	return final_adjusted_VAR



