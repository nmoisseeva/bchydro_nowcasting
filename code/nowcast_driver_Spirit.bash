#LOCAL BASH DRIVER FOR EMWXNET NOWCASTING SYSTEM ON SPIRIT
#==============================================================

#February 2016
#nmoisseeva@eos.ubc.ca

#--------------------------input-------------------------------
fcst_drive="/Volumes/forecasts/WRF3_RAW"
#-----------------------end of input---------------------------


echo "=============================================================="
echo "Initializing main nowcasting driver"


#get user defined variables from da_config.py
emx_dir=$(python - <<END
from da_config_spirit import *
print emx_dir 
END
)
netcdf_dir=$(python - <<END
from da_config_spirit import *
print netcdf_dir 
END
)
fig_dir=$(python - <<END
from da_config_spirit import *
print fig_dir
END
)
netcdf_prefix=$(python - <<END
from da_config_spirit import *
print netcdf_prefix
END
)
delay_hr=$(python - <<END
from da_config_spirit import *
print delay_hr
END
)
back_delay_hr=$(python - <<END
from da_config_spirit import *
print delay_hr + 1
END
)

#get current date and time in UTC, subtracting user-set time delay
# timestamp=$(date -u -v-${delay_hr}H +%Y%m%d%H)
timestamp=$(date --date="${delay_hr} hours ago" +%Y%m%d%H)
year=${timestamp:0:4}
yr=${timestamp:2:2}
month=${timestamp:4:2}
day=${timestamp:6:2}
hour=${timestamp:8:2}
fcst_init="00"
# old_stamp=$(date -u -v-${back_delay_hr}H +%Y%m%d%H)
old_stamp=$(date --date="${back_delay_hr} hours ago" +%Y%m%d%H)
old_hr=${old_stamp:8:2}
echo "Data assimilation will be performed for: $year-$month-$day $hour:00:00 "

#always pull new EmWxNet data (in case it is more complete)
emx_path=$emx_dir$year/$month/$day/	
mkdir -p $emx_path
echo "Extracting observations from the EmWxNet database"
../../emx_api/getstations_da.exe $year$month$day > /dev/null						
echo "Moving obs station data to $emx_path"
mv ./*.txt $emx_path > /dev/null	
 
# #check if the necessary NetCDF file aready exists
# netcdf_name=$netcdf_prefix$year-$month-$day"_"$hour:$fcst_init:00
# old_netcdf_name=$netcdf_prefix$year-$month-$day"_"$old_hr:$fcst_init:00
# netcdf_path="$netcdf_dir$year/$month/$day/$netcdf_name"
# old_netcdf_path="$netcdf_dir$year/$month/$day/$old_netcdf_name"
# if [ -e "$netcdf_path" ] && [ -e "$old_netcdf_path" ]; then
# 	echo "Found existing NetCDF files at: $netcdf_path"
# elif [ -e "$netcdf_path" ] && [ ! -e "$old_netcdf_path" ]; then
# 	echo "Downloading supporting NetCDF model data from Asiaq to: $old_netcdf_path"
# 	cp $fcst_drive/$yr$month$day$fcst_init/bsc-west/$old_netcdf_name $old_netcdf_path
# else
# 	echo "Moving required NetCDF model data from Asiaq to: $netcdf_path"
# 	mkdir -p $netcdf_dir$year/$month/$day/
# 	cp $fcst_drive/$yr$month$day$fcst_init/bsc-west/$netcdf_name $netcdf_path
# 	cp $fcst_drive/$yr$month$day$fcst_init/bsc-west/$old_netcdf_name $old_netcdf_path
# fi

# echo "Initializing main Python routine for data assimilation"
# mkdir -p $fig_dir/$year/$month/$day/$hour
# python main.py $netcdf_name

echo "=============================================================="
