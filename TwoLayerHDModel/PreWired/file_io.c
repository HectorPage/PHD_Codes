#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "file_io.h"

FILE* file_open(const char *fname, const char *format)
{
	FILE *fptr;
	
	fptr = fopen(fname, format);
	
	if(!fptr)
	{
		fprintf(stderr, "\nCould not open file %s\n", fname);
		exit(1);
	}
	
	return fptr;
}


void trim_space(char *str)
{
	char *end;
	
	end = str + (strlen(str)-1);
	
	while(isspace(*str))
		str++;
	
	while((end>=str)&&(isspace(*end)))
		end--;
	
	*(++end) = '\0';
	
	return;
}


void save_rates_activations(data *dptr, params *pptr)
{
	FILE *fptr;
	
	/* Save HD Rates */
	
	fptr = file_open("HDRates.bdat", "wb");
	
	fwrite(&(dptr->rates_HD_time[0][0]), sizeof(float), pptr->num_HD_cells*(pptr->timesteps/100), fptr);
	
	fclose(fptr);
	
	printf("HD Rates saved...\n");
	fflush(stdout);

	
	/* Save HD Activations */
	
	fptr = file_open("HDActivations.bdat", "wb");
	
	fwrite(&(dptr->activations_HD_time[0][0]), sizeof(float), pptr->num_HD_cells*(pptr->timesteps/100), fptr);
	
	fclose(fptr);
	
	printf("HD Activations saved...\n\n");
	fflush(stdout);

	
	/*Save COMB Rates*/
	
	fptr = file_open("COMBRates.bdat", "wb");
	
	fwrite(&(dptr->rates_COMB_time[0][0]), sizeof(float), pptr->num_COMB_cells*(pptr->timesteps/100), fptr);
	
	fclose(fptr);
	
	printf("COMB Rates saved...\n");
	fflush(stdout);

	
	/*Save COMB Activations*/
	
	fptr = file_open("COMBActivations.bdat", "wb");
	
	fwrite(&(dptr->activations_COMB_time[0][0]), sizeof(float), pptr->num_COMB_cells*(pptr->timesteps/100), fptr);
	
	fclose(fptr);
	
	printf("COMB Activations saved...\n");
	fflush(stdout);
	
	/*Save ROT Rates*/
	
	fptr = file_open("ROTRates.bdat", "wb");
	fwrite(&(dptr->rates_ROT_time[0]), sizeof(float), pptr->num_ROT_cells*(pptr->timesteps/100), fptr);
				
	fclose(fptr);
	
	printf("ROT Rates saved...\n");
	fflush(stdout);
	
	/*Save NOROT Rates*/
	
	fptr = file_open("NOROTRates.bdat", "wb");
	fwrite(&(dptr->rates_NOROT_time[0]), sizeof(float), pptr->num_NOROT_cells*(pptr->timesteps/100), fptr);
	
	fclose(fptr);
	
	printf("NOROT Rates saved...\n");
	fflush(stdout);

	
	return;
}

void save_pvector(data *dptr)
{
	FILE *fptr;
	
	fptr = file_open("Pvector.dat", "w");
	fprintf(fptr, "Conflict population vector is: %f ", dptr->pvector);
	
	fclose(fptr);
	
	return;
}


void save_initialiser(data *dptr, params *pptr)
{
	FILE *fptr;
	
	fptr = file_open("Input_Location1.bdat", "wb");
	fwrite(&(dptr->initialiser_location_time[0][0]), sizeof(float), pptr->num_HD_cells*(pptr->timesteps), fptr);

	
	fclose(fptr);
	
	return;
}


