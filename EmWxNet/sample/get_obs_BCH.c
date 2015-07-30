#include <stdio.h>
int main(int argc, char* argv[]) {

    int station=3477;  // call BC Hydro Edmonds station
    int st_date=20130101;
    int st_time=0;
    int end_date=20140322;
    int end_time=0;
    int level=0;

    // get_obs_ arrays for air temperature
    char varT[]="SFCTC";    // air temperature (degC)
    int* datesT;
    int* timesT;
    int sizeT = 0;
    float* obsT;
    // get_obs_ arrays for relative humidity
    char varRH[]="SFCRH";   // relative humidity (%)
    int* datesRH;
    int* timesRH;
    int sizeRH = 0;
    float* obsRH;
    // get_obs_ arrays for hourly precip 
    char varPPN[]="PCPTOT";     // hourly precip (mm)
    int* datesPPN;
    int* timesPPN;
    int sizePPN = 0;
    float* obsPPN;
    // get_obs_ arrays for wind speed
    char varWS[]="SFCWSPD";     // wind speed (km/h)
    int* datesWS;
    int* timesWS;
    int sizeWS = 0;
    float* obsWS;
    // get_obs_ arrays for wind direction 
    char varWD[]="SFCWDIR";     // wind direction (deg azimuth)
    int* datesWD;
    int* timesWD;
    int sizeWD = 0;
    float* obsWD;
    // get_obs_ arrays for wind chill
    char varWC[]="SFCWC";       // wind chill (degC)
    int* datesWC;
    int* timesWC;
    int sizeWC = 0;
    float* obsWC;
    // get_obs_ arrays for wind gust speed
    char varGS[]="SFCGS";       // wind gust speed (km/h)
    int* datesGS;
    int* timesGS;
    int sizeGS = 0;
    float* obsGS;
    // get_obs_ arrays for solar radiation 
    char varDSR[]="SFCSR";     // downwelling solar radiation (W/m^2)
    int* datesDSR;
    int* timesDSR;
    int sizeDSR = 0;
    float* obsDSR;
    // get_obs_ arrays for station pressure 
    char varP[]="SP";       // station pressure (hPa)
    int* datesP;
    int* timesP;
    int sizeP = 0;
    float* obsP;
    // get_obs_ arrays for maximum air temperature 
    char varTmax[]="SFCXTC";       // max air temp for archiving period (15 mins) (degC)
    int* datesTmax;
    int* timesTmax;
    int sizeTmax = 0;
    float* obsTmax;
    // get_obs_ arrays for minimum air temperature 
    char varTmin[]="SFCNTC";       // min air temp for archiving period (15 mins) (degC)
    int* datesTmin;
    int* timesTmin;
    int sizeTmin = 0;
    float* obsTmin;
    int k;
    
    // for output
    FILE* outputFile1a;    // Output ASCII file for UBC station
    char  outputFileName1a[100] = "BCHa.txt";    
    FILE* outputFile1b;    // Output ASCII file for UBC station
    char  outputFileName1b[100] = "BCHb.txt";    
    
    //****************************
    // Retrieve data for station*
    // ***************************
    sizeT = get_obs_(&station,&st_date,&st_time,&end_date,&end_time,&level,varT,&datesT,&timesT,&obsT,NULL,NULL);
    sizeRH = get_obs_(&station,&st_date,&st_time,&end_date,&end_time,&level,varRH,&datesRH,&timesRH,&obsRH,NULL,NULL);
    sizePPN = get_obs_(&station,&st_date,&st_time,&end_date,&end_time,&level,varPPN,&datesPPN,&timesPPN,&obsPPN,NULL,NULL);
    sizeWS = get_obs_(&station,&st_date,&st_time,&end_date,&end_time,&level,varWS,&datesWS,&timesWS,&obsWS,NULL,NULL);
    sizeWD = get_obs_(&station,&st_date,&st_time,&end_date,&end_time,&level,varWD,&datesWD,&timesWD,&obsWD,NULL,NULL);
    sizeWC = get_obs_(&station,&st_date,&st_time,&end_date,&end_time,&level,varWC,&datesWC,&timesWC,&obsWC,NULL,NULL);
    sizeGS = get_obs_(&station,&st_date,&st_time,&end_date,&end_time,&level,varGS,&datesGS,&timesGS,&obsGS,NULL,NULL);
    sizeDSR = get_obs_(&station,&st_date,&st_time,&end_date,&end_time,&level,varDSR,&datesDSR,&timesDSR,&obsDSR,NULL,NULL);
    sizeP = get_obs_(&station,&st_date,&st_time,&end_date,&end_time,&level,varP,&datesP,&timesP,&obsP,NULL,NULL);
    sizeTmax = get_obs_(&station,&st_date,&st_time,&end_date,&end_time,&level,varTmax,&datesTmax,&timesTmax,&obsTmax,NULL,NULL);
    sizeTmin = get_obs_(&station,&st_date,&st_time,&end_date,&end_time,&level,varTmin,&datesTmin,&timesTmin,&obsTmin,NULL,NULL);

    //*************************
    // Open files to write to *
    //*************************
    outputFile1a = fopen(outputFileName1a, "w+");   // station 1: air temp, RH, ppn, wind speed, gust speed
    fprintf(outputFile1a, " %s %s %s %s %s %s %s\n", "Dates","Times","Temp (degC)","RH (%)","Hourly precip (mm)","Wind speed (km/h)","Gust speed (km/h)");
    
    outputFile1b = fopen(outputFileName1b, "w+");   // station 1: wind direction, wind chill, solar radiation, pressure
    fprintf(outputFile1b, " %s %s %s %s %s %s %s %s\n", "Dates","Times","Wind dir (deg)","Wind chill (degC)","Downwelling solar radiation (W/m^2)","Station pressure (hPa)","Tmax (degC)","Tmin (degC)");
    if(!outputFile1a || !outputFile1b) 
    {
        // Could not open file
        printf("get_allobs.c: WARNING: could not open output file\n");
    }
    
    // Write to file
    if(sizeT == 0)
    {
        printf("    get_allobs.c: Warning: No observations retrieved\n");
        printf("    sizeT = %i\n", sizeT);
    }
    else {
        for( k=0; k<sizeT; k++ )
        {   
            //printf("    sizeT1 = %i\n", sizeT1);
            fprintf(outputFile1a, " %d %d %6.2f %6.2f %6.2f %6.2f %6.2f\n", datesT[k], timesT[k], obsT[k], obsRH[k], obsPPN[k], obsWS[k], obsGS[k]);
            fprintf(outputFile1b, "  %d %d %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n", datesT[k], timesT[k], obsWD[k], obsWC[k], obsDSR[k], obsP[k], obsTmax[k], obsTmin[k]);
        }
    }

    fclose(outputFile1a);
    fclose(outputFile1b);

    // clean up
    free(datesT);
    free(datesRH);
    free(datesPPN);
    free(datesWS);
    free(datesWD);
    free(datesWC);
    free(datesGS);
    free(datesDSR);
    free(datesP);
    free(datesTmax);
    free(datesTmin);
    
    free(timesT);
    free(timesRH);
    free(timesPPN);
    free(timesWS);
    free(timesWD);
    free(timesWC);
    free(timesGS);
    free(timesDSR);
    free(timesP);
    free(timesTmax);
    free(timesTmin);
    
    free(obsT);
    free(obsRH);
    free(obsPPN);
    free(obsWS);
    free(obsWD);
    free(obsWC);
    free(obsGS);
    free(obsDSR);
    free(obsP);
    free(obsTmax);
    free(obsTmin);
    
    return 0;
}
