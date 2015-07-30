# CONFIGURATION FILE FOR EMWXNET NOWCASTING DATA ASSIMILATION

# nmoisseeva@eos.ubc.ca
# June 2015


#------------------------------------------------------
#define input and output locations and source files
fig_dir = '/Users/nadya2/Workspace/nowcasting/figs/'		#directory for storing figures
netcdf_dir = '/Users/nadya2/data/netcdf/'					#directory of raw netcdf data
emx_dir = '/Users/nadya2/data/emwxnet/20150224/'			#directory of EmWeatherNet data
emx_name = 'selectStnList.txt'								#output file from getstsations_da.c 
geotiff = '12arcsecDEM.tif'									#DEM model geoTiff file
landmask_dir = '/Users/nadya2/Workspace/nowcasting/npy/'	#folder containing landmask (optional)
basemap_dir = '/Users/nadya2/Workspace/nowcasting/npy/'		#folder containing basemap file (optional)

netcdf_name = 'wrfout_d03_2015-02-24_12:00:00'				#name of raw netcdf file

#------------------------------------------------------
#domain configuration
lat = [48.2,51.]				#[min lat, max lat]
lon = [-128.8,-121.]			#[min lon, max lon]

#lapse rate calculation adjustments
lvl = 20 						#max model level to use for lapse rate calculation 

#------------------------------------------------------
#data assimilation configuration
run_MD = 1 						#flag to run mother-daughter (alternatively ROI-based)
verif_frac = 0.7				#fraction of data to use for training (remainder for verification)
bias_mode = 0. 					#flag to correct temperature increments (alternatively will correct temperature)



#------------------------------------------------------
#ROI DA configuration
roi = 0.09						#horizontal influence range (decimal degrees)
elev_roi = 500. 				#vertical influence range (meters)

#------------------------------------------------------
#MD DA configuration
params = 2,2,750,500			#[a,b, Z_ref1, Z_ref2] - parameters for sharing factor calculation
dist_cutoff = 1. 				#maximum anisotropic horizontal distance in degrees to continue iteration
weight_cutoff = 0.01			#minimum change in weight values required to continue iteration



#ADDITIONS TO TEST:
#correct bias (increment) - no actual temperature
#correct colorbar height
#make flexible colormap arange
#test using MD distance with inverse