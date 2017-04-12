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
	
	fptr = file_open("HDRates.bdat", "wb");
	
	fwrite(&(dptr->rates_E1_time[0][0]), sizeof(float), pptr->num_E1_cells*(pptr->timesteps/100), fptr);
	
	fclose(fptr);
	

	/* Save E1 Activations */
	
	fptr = file_open("HDActivations.bdat", "wb");
	
	fwrite(&(dptr->activations_E1_time[0][0]), sizeof(float), pptr->num_E1_cells*(pptr->timesteps/100), fptr);
	
	fclose(fptr);
	
		
	return;
}


void save_input(data *dptr, params *pptr)
{
	FILE *fptr;
	
	fptr = file_open("Input_Location1.bdat", "wb");
	fwrite(&(dptr->input_location_time[0][0]), sizeof(float), pptr->num_E1_cells*(pptr->timesteps/100), fptr);

	
	fclose(fptr);
	
	return;
}


void load_full_RC(connect *cptr, params *pptr)
{
	int post,pre;
	int synapse;
	
	for(post=0; post<pptr->num_E1_cells; post++)
	{
	
		for(pre=0; pre<pptr->num_E1_cells; pre++)
		{
			
			cptr->full_RC[post][pre] = 0.0;
			
		}
	}
	
		
	for(post=0; post<pptr->num_E1_cells; post++)
	{
		for(pre=0; pre<cptr->num_RC1_connections; pre++)
		{
			synapse = cptr->RC1_connections[post][pre];
			
			cptr->full_RC[post][synapse] = cptr->RC1_weights[post][pre];
			
		}
	}
	
}


void save_RCweights_init(connect *cptr, params *pptr)
{
	FILE *fptr;
	
	fptr = file_open("initial_RCweights.bdat","wb");
	fwrite(&(cptr->full_RC[0][0]), sizeof(float), pptr->num_E1_cells*pptr->num_E1_cells, fptr);
	
	fclose(fptr);
	
	return;
	
}


void save_RCweights_intermediate(connect *cptr, params *pptr)
{
	FILE *fptr;
	
	fptr = file_open("middle_RCweights.bdat","wb");
	fwrite(&(cptr->full_RC[0][0]), sizeof(float), pptr->num_E1_cells*pptr->num_E1_cells, fptr);
	
	fclose(fptr);
	
	return;
	
}

void save_RCweights_final(connect *cptr, params *pptr)
{
	FILE *fptr;
	
	fptr = file_open("final_RCweights.bdat","wb");
	fwrite(&(cptr->full_RC[0][0]), sizeof(float), pptr->num_E1_cells*pptr->num_E1_cells, fptr);
	
	fclose(fptr);
	
	return;
	
}

void save_weight_vectors(data *dptr, params *pptr)
{
	FILE *fptr;

	fptr = file_open("intermediate_weight_vectors.bdat", "wb");
	fwrite(&(dptr->weight_vector_intermediate[0]),sizeof(float), pptr->num_E1_cells, fptr);
	
	fclose(fptr);
	
	fptr = file_open("final_weight_vectors.bdat", "wb");
	fwrite(&(dptr->weight_vector_final[0]),sizeof(float), pptr->num_E1_cells, fptr);
	
	fclose(fptr);
	
	return;
	
}

void save_wVecs(data *dptr, params *pptr)
{
	FILE *fptr;
	
	fptr = file_open("weight_vectors.bdat", "wb");
	fwrite(&(dptr->wVec_time[0][0]), sizeof(float), pptr->num_E1_cells*(pptr->timesteps/100), fptr);
	
	fclose(fptr);
	
	return;
	
}


