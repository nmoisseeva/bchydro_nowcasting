# CONFIGURATION FILE FOR EMWXNET NOWCASTING DATA ASSIMILATION

# nmoisseeva@eos.ubc.ca
# June 2015


#------------------------------------------------------
#defining input and output locations and source files
fig_dir = './figs/'											#directory for storing figures
netcdf_dir = '/Users/nadya2/data/netcdf/'					#directory of raw netcdf data
netcdf_prefix = 'wrfout_d03_'								#prefix format of raw NetCDF file
emx_dir = '/Users/nadya2/data/emwxnet/'						#directory of EmWeatherNet data
emx_name = 'selectStnList.txt'								#output file from getstsations_da.c 
geotiff = '12arcsecDEM.tif'									#DEM model geoTiff file

landmask_dir = './npy/'										#folder containing landmask (optional, no changes recommended)
basemap_dir = './npy/'										#folder containing basemap file (optional, no changes recommended)

#------------------------------------------------------
#domain configuration
lat = [48.2,51.]				#[min lat, max lat]
lon = [-128.8,-121.]			#[min lon, max lon]

#------------------------------------------------------
#data assimilation configuration

run_MD = 1 						#flag to run mother-daughter (alternatively ROI-based)
verif_frac = 0.7				#fraction of data to use for training (remainder for verification)
bias_mode = 0 					#flag to correct temperature increments (alternatively will correct temperature)
delay_hr = 4 	 				#for operational use: analysis delay (in hours)

#------------------------------------------------------
#ROI DA configuration
roi = 0.09						#horizontal influence range (decimal degrees)
elev_roi = 500. 				#vertical influence range (meters)

#------------------------------------------------------
#MD DA configuration
params = 2,2,750,650			#a,b, Z_ref1, Z_ref2 - parameters for sharing factor calculation
dist_cutoff = 1. 				#maximum anisotropic horizontal distance in degrees to continue iteration
weight_cutoff = 0.01			#minimum change in weight values required to continue iteration

#------------------------------------------------------
#temperature and lapse rate configuration
lvl = 20 						#max model level to use for lapse rate calculation 
T_range = [-50,50]				#colormap range for Temperature plots in C (center at 0)



#BUGS to fix:
#fix keychain issues on Spirit
#flexible graphics output size
#md - save max dist tag
#md _ md weights-dist plot - where to save?
#colorcode stations temperatures
#bc hydro tags