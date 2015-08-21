
//====================== block_pull.c ===================
//nmoisseeva@eos.ubc.ca
//February 2015

//script pulls station info for defined lon/lat box and select date range,
//then stores their data for select vars as text files

//=======================================================
#include <stdio.h>
#include <stdlib.h>

 int main(int argc, char *argv[] ) 
 {
	//Part (1)
	//pull all station id's and info for defined gridbox
	//----------------------- iput (1) ------------------------
	// int st_date = 20150224;						//start date
	// int end_date = 20150225;					//end date
	int st_date = atoi(argv[1]);						//start date
	int end_date = st_date+1;					//end date
	float min_lat = 48.0000;					//define box boundaries
	float max_lat = 52.0000;
	float min_long = -129.0000;
	float max_long = -120.0000;
	int size = 2000;							//max array size (for all EmWxNet ~4000)
	char outputID[] = "selectStnList.txt";		//output file

	//------------------- end of input (1) --------------------

	printf("\nPart (1): Extract all station within the defined box domain:\n");

	//decare variables and pointers used by API
	int numStat, nStn;
	int* station_ids, *elevation;	
	float* latitude, *longitude;
	FILE* stnListoutput;

	//allocate space and run API script
	station_ids = (int*)malloc(sizeof(int) * size);
	elevation = (int*)malloc(sizeof(int) * size);
	latitude = (float*)malloc(sizeof(float) * size);
	longitude = (float*)malloc(sizeof(float) * size);
	numStat = getstations_(&st_date, &end_date, &min_lat, &max_lat, &min_long, &max_long, station_ids, latitude, longitude, elevation, &size, NULL, NULL);


	//test data and write it to file
	stnListoutput = fopen(outputID, "w+");				
	fprintf(stnListoutput, "%s %s %s %s \n", "Station ID", "Latitude", "Longitude", "Elevation");
	if(!stnListoutput) 									
	{
		printf("ERROR: could not open output file %s \n", outputID);
	}
	if(numStat == 0)									
	{
		printf("WARNING: no stations matched the selected criteria\n");
	}
	else 
	{
		for(nStn=0; nStn<numStat; nStn++)				
		{
			fprintf(stnListoutput, "%d %6.2f %6.2f %d \n" , station_ids[nStn],latitude[nStn],longitude[nStn],elevation[nStn]);
		}
	}
	fclose(stnListoutput);
	printf("%d stations extracted from EmWxNet\n", numStat);
	printf("------------------ Part (1) complete---------------\n");



	//Part (2)
	//get observations for the selected station_ids
	//----------------------- input (2) ------------------------
	char varT[] = "SFCTC";
	char varWS[] = "SFCWSPD";
	char varWD[] = "SFCWDIR";
	float* valuesT, *valuesWS, *valuesWD;
	int numT, numWS, numWD;
	int* datesT, *datesWS, *datesWD;
	int* timesT, *timesWS, *timesWD;


	//------------------- end of input (2) --------------------

	printf("\nPart (2): Extract observations for the retrieved stations:\n");

	//decare variables and pointers used by API
	int st_time = 0;
	int end_time = 0;
	int level = 0;
	int station, nObs;  
	FILE* obsFile;

	for(nStn=0; nStn<numStat; nStn++)
	{
		//create individual text file with data for each station
		station = station_ids[nStn];
		printf("Station: %d \n", station);
		char obsFilename[100];
		snprintf(obsFilename, sizeof(char)*100, "%d.txt",station);

		numT = get_obs_(&station,&st_date,&st_time,&end_date,&end_time,&level,varT,&datesT,&timesT,&valuesT,NULL,NULL);
		numWS = get_obs_(&station,&st_date,&st_time,&end_date,&end_time,&level,varWS,&datesWS,&timesWS,&valuesWS,NULL,NULL);
		numWD = get_obs_(&station,&st_date,&st_time,&end_date,&end_time,&level,varWD,&datesWD,&timesWD,&valuesWD,NULL,NULL);


		if(numT == 0 && numWS == 0 && numWD== 0)									
			{
				printf("WARNING: no observations matched the selected criteria\n");
			}
	
		else 
			{
				obsFile = fopen(obsFilename,"w+");
				fprintf(obsFile,"%s %s %s %s %s \n", "Date", "Time", "Temperature (degC)","Wind Speed (km/h)","Wind Dir (deg)");
				if(!obsFile) 									
				{
					printf("ERROR: could not open output file %s \n", obsFilename);
				}
				for(nObs=0; nObs<numT; nObs++)				
				{
					fprintf(obsFile, "%d %d %6.2f %6.2f %6.2f \n" ,datesT[nObs],timesT[nObs],valuesT[nObs],valuesWS[nObs],valuesWD[nObs]);
				}
				fclose(obsFile);

			}

	}

	printf("Data extracted from %d EmWxNet stations\n", numStat);
	printf("------------------ Part (2) complete---------------\n");


	free(station_ids);
	free(elevation);
	free(latitude);
	free(longitude);

 }


