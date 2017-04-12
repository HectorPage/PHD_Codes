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
	
	/* Save E1 Rates */
	
	fptr = file_open("E1Rates.bdat", "wb");
	
	fwrite(&(dptr->rates_E1_time[0][0]), sizeof(float), pptr->num_E1_cells*(pptr->timesteps/100), fptr);
	
	fclose(fptr);
	

	/* Save E1 Activations */
	
	fptr = file_open("E1Activations.bdat", "wb");
	
	fwrite(&(dptr->activations_E1_time[0][0]), sizeof(float), pptr->num_E1_cells*(pptr->timesteps/100), fptr);
	
	fclose(fptr);
	
		
	return;
}

void save_pvector(data *dptr)
{
	FILE *fptr;
	
	fptr = file_open("Pvector.dat", "w");
	fprintf(fptr, "End population vector is: %f ", dptr->pvector);
	
	fclose(fptr);
	
	return;
}

void save_visinput1(data *dptr, params *pptr)
{
	FILE *fptr;
	
	fptr = file_open("Input_Location1.bdat", "wb");
	fwrite(&(dptr->vis1_location_time[0][0]), sizeof(float), pptr->num_E1_cells*(pptr->timesteps/100), fptr);

	
	fclose(fptr);
	
	return;
}

void save_visring(data *dptr, params *pptr)
{
	FILE *fptr;
	
	fptr = file_open("visring_location.bdat","wb");
	fwrite(&(dptr->visring_location_time[0][0]), sizeof(float), pptr->num_E1_cells*(pptr->timesteps/100), fptr);
	
	fclose(fptr);
	
	fptr = file_open("i_total.bdat","wb");
	fwrite(&(dptr->i_total[0]), sizeof(float), pptr->num_E1_cells, fptr);
	
	fclose(fptr);	
	
	return;
	
}


void save_ffweights_init(connect *cptr, params *pptr)
{
	FILE *fptr;
	
	fptr = file_open("initial_ffweights.bdat","wb");
	fwrite(&(cptr->weights_visring[0][0]), sizeof(float), pptr->num_E1_cells*pptr->num_E1_cells, fptr);
	
	fclose(fptr);
	
	return;
	
}


void save_ffweights_final(connect *cptr, params *pptr)
{
	FILE *fptr;
	
	fptr = file_open("final_ffweights.bdat","wb");
	fwrite(&(cptr->weights_visring[0][0]), sizeof(float), pptr->num_E1_cells*pptr->num_E1_cells, fptr);
	
	fclose(fptr);
	
	return;
	
}

void weight_dump(connect *cptr, params *pptr)
{
	FILE *fptr;
	
	fptr = file_open("weight_dump.bdat","ab");
	fwrite(&(cptr->weights_visring[0][0]), sizeof(float), pptr->num_E1_cells*pptr->num_E1_cells, fptr);
	
	fclose(fptr);
	
	return;
	
}
