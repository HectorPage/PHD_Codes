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
	
	printf("\n\nHD Rates saved...\n\n");
	fflush(stdout);

	
	/* Save HD Activations */
	
	fptr = file_open("HDActivations.bdat", "wb");
	
	fwrite(&(dptr->activations_HD_time[0][0]), sizeof(float), pptr->num_HD_cells*(pptr->timesteps/100), fptr);
	
	fclose(fptr);
	
	printf("\n\nHD Activations saved...\n\n");
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

void save_weight_vector(data *dptr, params *pptr)
{
	FILE *fptr;
	
	fptr = file_open("initial_ff_weight_vectors.bdat", "wb");
	fwrite(&(dptr->HD_weight_vectors_init[0]), sizeof(float), pptr->num_HD_cells, fptr);
	fclose(fptr);
	
	
	fptr = file_open("intermediate_ff_weight_vectors.bdat", "wb");
	fwrite(&(dptr->HD_weight_vectors_intermediate[0]), sizeof(float), pptr->num_HD_cells, fptr);
	fclose(fptr);
	
	
	fptr = file_open("final_ff_weight_vectors.bdat", "wb");
	fwrite(&(dptr->HD_weight_vectors_final[0]), sizeof(float), pptr->num_HD_cells, fptr);
	fclose(fptr);
	
	return;
	
}

void save_pi_input(data *dptr, params *pptr)
{
	FILE *fptr;
	
	fptr = file_open("Input_Location1.bdat", "wb");
	fwrite(&(dptr->pi_location_time[0][0]), sizeof(float), pptr->num_HD_cells*(pptr->timesteps/100), fptr);

	
	fclose(fptr);
	
	return;
}

void save_visring(data *dptr, params *pptr)
{
	FILE *fptr;
	
	fptr = file_open("visring_location.bdat","wb");
	fwrite(&(dptr->visring_location_time[0][0]), sizeof(float), pptr->num_HD_cells*(pptr->timesteps/100), fptr);
	
	fclose(fptr);
		
	return;
	
}


void save_ffweights_init(connect *cptr, params *pptr)
{
	FILE *fptr;
	
	fptr = file_open("initial_ffweights.bdat","wb");
	fwrite(&(cptr->weights_visring[0][0]), sizeof(float), pptr->num_HD_cells*pptr->num_HD_cells, fptr);
	
	fclose(fptr);
	
	return;
	
}


void save_ffweights_intermediate(connect *cptr, params *pptr)
{
	FILE *fptr;
	
	fptr = file_open("intermediate_ffweights.bdat","wb");
	fwrite(&(cptr->weights_visring[0][0]), sizeof(float), pptr->num_HD_cells*pptr->num_HD_cells, fptr);
	
	fclose(fptr);
	
	return;
	
}

void save_ffweights_final(connect *cptr, params *pptr)
{
	FILE *fptr;
	
	fptr = file_open("final_ffweights.bdat","wb");
	fwrite(&(cptr->weights_visring[0][0]), sizeof(float), pptr->num_HD_cells*pptr->num_HD_cells, fptr);
	
	fclose(fptr);
	
	return;
	
}


void save_weight_vectors(data *dptr, params *pptr)
{
	FILE *fptr;
	
	fptr = file_open("initial_weight_vectors.bdat", "wb");
	fwrite(&(dptr->HD_weight_vectors_init[0]),sizeof(float), pptr->num_HD_cells, fptr);
	
	fclose(fptr);
	
	fptr = file_open("intermediate_weight_vectors.bdat", "wb");
	fwrite(&(dptr->HD_weight_vectors_intermediate[0]),sizeof(float), pptr->num_HD_cells, fptr);
	
	fclose(fptr);
	
	fptr = file_open("final_weight_vectors.bdat", "wb");
	fwrite(&(dptr->HD_weight_vectors_final[0]),sizeof(float), pptr->num_HD_cells, fptr);
	
	fclose(fptr);
	
	
}
