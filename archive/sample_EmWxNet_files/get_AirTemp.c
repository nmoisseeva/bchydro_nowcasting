#include <stdio.h>
int main(int argc, char* argv[]) {

    int station=417;
    int st_date=20130101;
    int st_time=0;
    int end_date=20130701;
    int end_time=0;
    int level=0;
    char type[]="SFCTC";
    int* dates;
    int* times;
    float* values;
    int k;
    
    // for output
    FILE* outputFile;    // Output ASCII file
    char  outputFileName[100] = "2013AirTemp:417.txt";    
    // char  outputFileName[100] = "/users/model/src/rosie/solar_rad.txt";    

    int i = get_obs_(&station,&st_date,&st_time,&end_date,&end_time,&level,type,&dates,&times,&values,NULL,NULL);
    printf("we found %d obs\n",i);

    //for(k=0; k<i; k++)
    //{
    //    printf("  %d %d %f\n", dates[k], times[k], values[k] );
    //}

    // Open file to write to
    outputFile = fopen(outputFileName, "w+");
    if(!outputFile) {
        // Could not open file
        printf("get_AirTemp.c: WARNING: could not open output file\n");
    }
    
    for( k=0; k<i; k++ )
    {
        fprintf(outputFile, "  %d %d %6.2f\n", dates[k], times[k], values[k] );
    }

    fclose(outputFile);

    free(dates);
    free(times);
    free(values);

    return 0;
}
