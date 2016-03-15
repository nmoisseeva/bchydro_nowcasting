Contact info: Nadya Moisseeva (nmoisseeva@eos.ubc.ca)

-------------------------------------------------------------
PLEASE REFER TO ./docs/v2.0/UserManual.pdf FOR INSTRUCTIONS
-------------------------------------------------------------

***Description:***
WRF/EmWxNet High Resolution Nowcasting System produces refined analysis fields by combining model forecast with real-time observations through downscaling and data assimilation.The current system relies on meteorological fields from archived WRF runs and observations from Emergency Weather Net Database (EmWxNet). Graphics are saved in /nfs/neltharion/www/results/HRSA/YYMMDDHH/PNG/g6/

***Setup Requirements:***
Python 2.7 (including numpy, scipy, gdal, basemaps/mpl_toolkits, matplotlib)

***Operational Setup Instructions for Schig:***
1. Compile contents of ./emx subfolder with the provided makefile.
2. Edit the header of nowcasting_driver.bash in the main directory to point to the compiled getstations_da.exe
3. cron a bash call to nowcasting_driver.bash (e.g. 0 * * * * bash $PATH/nowcasting_driver.bash >> $PATH/`date +\%Y\%m\%d\%H\%M`-hrsa.log 2>&1)

***The above instructions are specific to operational setup: for all other usage refer to User Manual***
