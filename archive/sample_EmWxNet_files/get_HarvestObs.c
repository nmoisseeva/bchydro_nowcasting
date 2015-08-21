#include <stdio.h>
#include <stdlib.h>
int main(int argc, char* argv[]) {

    int station1=597;  // YVR
    int station2=396;  // North Delta
    int station3=400;  // Burnaby South
    int station4=399;  // Richmond South
    int st_date=20140801;
    int st_time=0;
    int end_date=20140901;
    int end_time=0;
    int level=0;

    // get_obs_ arrays for wind speed
    char varWS[]="SFCWSPD";     // wind speed (km/h)
    int* datesWS1;
    int* timesWS1;
    int sizeWS1 = 0;
    float* obsWS1;
    int* datesWS2;
    int* timesWS2;
    int sizeWS2 = 0;
    float* obsWS2;
    int* datesWS3;
    int* timesWS3;
    int sizeWS3 = 0;
    float* obsWS3;
    int* datesWS4;
    int* timesWS4;
    int sizeWS4 = 0;
    float* obsWS4;
    // get_obs_ arrays for wind direction 
    char varWD[]="SFCWDIR";     // wind direction (deg azimuth)
    int* datesWD1;
    int* timesWD1;
    int sizeWD1 = 0;
    float* obsWD1;
    int* datesWD2;
    int* timesWD2;
    int sizeWD2 = 0;
    float* obsWD2;
    int* datesWD3;
    int* timesWD3;
    int sizeWD3 = 0;
    float* obsWD3;
    int* datesWD4;
    int* timesWD4;
    int sizeWD4 = 0;
    float* obsWD4;
    
    int k;
    
    // for output
    FILE* outputFile1;    // Output ASCII file
    char  outputFileName1[100] = "Obs_Aug597.txt";    
    FILE* outputFile2;    // Output ASCII file
    char  outputFileName2[100] = "Obs_Aug396.txt";    
    FILE* outputFile3;    // Output ASCII file
    char  outputFileName3[100] = "Obs_Aug400.txt";    
    FILE* outputFile4;    // Output ASCII file
    char  outputFileName4[100] = "Obs_Aug399.txt";    
    
    //****************************
    // Retrieve data for station*
    // ***************************
    sizeWS1 = get_obs_(&station1,&st_date,&st_time,&end_date,&end_time,&level,varWS,&datesWS1,&timesWS1,&obsWS1,NULL,NULL);
    sizeWS2 = get_obs_(&station2,&st_date,&st_time,&end_date,&end_time,&level,varWS,&datesWS2,&timesWS2,&obsWS2,NULL,NULL);
    sizeWS3 = get_obs_(&station3,&st_date,&st_time,&end_date,&end_time,&level,varWS,&datesWS3,&timesWS3,&obsWS3,NULL,NULL);
    sizeWS4 = get_obs_(&station4,&st_date,&st_time,&end_date,&end_time,&level,varWS,&datesWS4,&timesWS4,&obsWS4,NULL,NULL);
    sizeWD1 = get_obs_(&station1,&st_date,&st_time,&end_date,&end_time,&level,varWD,&datesWD1,&timesWD1,&obsWD1,NULL,NULL);
    sizeWD2 = get_obs_(&station2,&st_date,&st_time,&end_date,&end_time,&level,varWD,&datesWD2,&timesWD2,&obsWD2,NULL,NULL);
    sizeWD3 = get_obs_(&station3,&st_date,&st_time,&end_date,&end_time,&level,varWD,&datesWD3,&timesWD3,&obsWD3,NULL,NULL);
    sizeWD4 = get_obs_(&station4,&st_date,&st_time,&end_date,&end_time,&level,varWD,&datesWD4,&timesWD4,&obsWD4,NULL,NULL);

    //*************************
    // Open files to write to *
    //*************************
    outputFile1 = fopen(outputFileName1, "w+");   // station 1: wind speed, wind direction
    outputFile2 = fopen(outputFileName2, "w+");   // station 2: wind speed, wind direction
    outputFile3 = fopen(outputFileName3, "w+");   // station 3: wind speed, wind direction
    outputFile4 = fopen(outputFileName4, "w+");   // station 4: wind speed, wind direction 
    if(!outputFile1 || !outputFile2 || !outputFile3 || !outputFile4) 
    {
        // Could not open file
        printf("get_HarvestObs.c: WARNING: could not open output file\n");
    }
    
    // Write to file
    if(sizeWS1 == 0 || sizeWS2 == 0 || sizeWS3 == 0 || sizeWS4 == 0 || sizeWD1 == 0 || sizeWD2 == 0 || sizeWD3 == 0 || sizeWD4 == 0)
    {
        printf("    get_HarvestObs.c: Warning: No observations retrieved\n");
        printf("    sizeWS1 = %i\n", sizeWS1);
        printf("    sizeWS2 = %i\n", sizeWS2);
        printf("    sizeWS3 = %i\n", sizeWS3);
        printf("    sizeWS4 = %i\n", sizeWS4);
        printf("    sizeWD1 = %i\n", sizeWD1);
        printf("    sizeWD2 = %i\n", sizeWD2);
        printf("    sizeWD3 = %i\n", sizeWD3);
        printf("    sizeWD4 = %i\n", sizeWD4);
    }
    else {
        for( k=0; k<sizeWS1; k++ )
        {   
            //printf("    sizeT1 = %i\n", sizeT1);
            fprintf(outputFile1, " %d %d %6.2f %6.2f\n", datesWS1[k], timesWS1[k], obsWS1[k], obsWD1[k]);
            fprintf(outputFile2, " %d %d %6.2f %6.2f\n", datesWS2[k], timesWS2[k], obsWS2[k], obsWD2[k]);
            fprintf(outputFile3, " %d %d %6.2f %6.2f\n", datesWS3[k], timesWS3[k], obsWS3[k], obsWD3[k]);
            fprintf(outputFile4, " %d %d %6.2f %6.2f\n", datesWS4[k], timesWS4[k], obsWS4[k], obsWD4[k]);
        }
    }

    fclose(outputFile1);
    fclose(outputFile2);
    fclose(outputFile3);
    fclose(outputFile4);

    // clean up
    free(datesWS1);
    free(datesWD1);
    free(datesWS2);
    free(datesWD2);
    free(datesWS3);
    free(datesWD3);
    free(datesWS4);
    free(datesWD4);
    
    free(timesWS1);
    free(timesWD1);
    free(timesWS2);
    free(timesWD2);
    free(timesWS3);
    free(timesWD3);
    free(timesWS4);
    free(timesWD4);
    
    free(obsWS1);
    free(obsWD1);
    free(obsWS2);
    free(obsWD2);
    free(obsWS3);
    free(obsWD3);
    free(obsWS4);
    free(obsWD4);
    
    return 0;
}
