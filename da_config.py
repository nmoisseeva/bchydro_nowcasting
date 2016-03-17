# CONFIGURATION FILE FOR EMWXNET NOWCASTING DATA ASSIMILATION

# nmoisseeva@eos.ubc.ca
# February 2016


#------------------------------------------------------
#defining input and output locations and source files
fig_dir = '/tmp/WAN00WP03/data/bc_nowcast/figs/'						#directory for storing figures
netcdf_dir = '/tmp/WAN00WP03/data/netcdf/'							#directory of raw netcdf data
netcdf_prefix = 'wrfout_d03_'								#prefix format of raw NetCDF file
emx_dir = '/tmp/WAN00WP03/data/emwxnet/'								#directory of EmWeatherNet data
emx_name = 'selectStnList.txt'								#output file from getstsations_da.c 
elev_geotiff = '12arcsecDEM.tif'							#DEM model geoTiff file
aspect_geotiff = '12arcsecASPECT.tif'						#aspect geoTiff file
slope_geotiff = '12arcsecSLOPE.tif'							#slope geoTiff file
landmask_dir = './npy/'										#folder containing landmask (optional, no changes recommended)
basemap_dir = './npy/'										#folder containing basemap file (optional, no changes recommended)


#------------------------------------------------------
#domain configuration
lat = [48.2,51.]				#[min lat, max lat]
lon = [-128.8,-121.]			#[min lon, max lon]


#------------------------------------------------------
#general settings
delay_hr = 3	 	 			#for operational use: analysis delay (in hours)
verif_frac = 0.7				#fraction of data to use for training (remainder for verification)


#------------------------------------------------------
#temperature assimilation configuration
temp_run_MD = 1 				#flag to run mother-daughter for temperature (alternatively ROI-based)
temp_bias_mode = 1 				#flag to correct temperature increments (alternatively will correct temperature)
temp_roi = 0.2					#horizontal influence range (decimal degrees) for temperature
temp_elev_roi = 500. 			#vertical influence range for temperature(meters)
lvl = 20 						#max model level to use for lapse rate calculation 
T_range = [-35,35]				#colormap range for Temperature plots in C (center at 0)
#MD DA configuration
params = 2,2,750,650			#a,b, Z_ref1, Z_ref2 - parameters for sharing factor calculation
dist_cutoff = 1.0 				#maximum anisotropic horizontal distance in degrees to continue iteration (degrees to 1 decimal)
weight_cutoff = 0.01			#minimum change in weight values required to continue iteration


#------------------------------------------------------
#precipitation assimilation configuration
precip_run_PRISM = 1 			#flag to run prism assimilation for precip (alternatively ROI-based)
precip_bias_mode = 1			#flag to correct temperature increments (alternatively will correct temperature)
rain_roi = 0.2 	 	 			#horizontal influence range (decimal degrees) for rain
precip_elev_roi = 2000. 		#vertical influence range for precipitations (meters)
R_range = [0.5,100] 				#colormap range for rain (will be plotted on log scale)


# #------------------------------------------------------
# #wind assimilation configuration
# wind_run_MD = 0 				#flag to run mother-daughter for wind (ROI-based is recommended!)
# wind_bias_mode = 0 				#flag to correct wind increments (alternatively will correct temperature)
# wind_roi = 0.2 					#horizontal influence range (decimal degrees) for wind
# wind_elev_roi = 1000. 			#vertical influence range for wind (meters)
# W_range = [-20,20]  			#colormap range for u- and v- wind plots (center at 0)









#BUGS to fix:
#change colormap for prism weght plots to start on white
#make glossary in manual
#remove landmask calls from MD.py (leave unmasked)
