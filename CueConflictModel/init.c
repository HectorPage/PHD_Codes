#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_randist.h>

#include "init.h"
#include "data.h"
#include "file_io.h"

void read_parameters(const char *fname, params *pptr, connect *cptr)
{
	FILE *fptr;
	char str[200], *tokptr;
	
	fptr = file_open(fname, "r");
	
	while(!feof(fptr))
	{
		fgets(str, 200, fptr);
		
		tokptr = strtok(str, "=");
		
		while(tokptr!=NULL)
		{
			trim_space(tokptr);
			
			if(strcmp(tokptr, "time")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->time = atof(tokptr);
			}
			else if(strcmp(tokptr, "E1 cells")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->num_E1_cells = atoi(tokptr);
			}
			else if(strcmp(tokptr, "RC1 connections")==0)
			{
				tokptr = strtok(NULL, "=");
				cptr->num_RC1_connections = atoi(tokptr);
			}
			else if(strcmp(tokptr, "phi RC1")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->phi_RC1 = atof(tokptr);
			}
			else if(strcmp(tokptr, "tau E1")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->tau_E1 = atof(tokptr);
			}
			else if(strcmp(tokptr, "alpha E1")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->alpha_E1 = atof(tokptr);
			}
			else if(strcmp(tokptr, "beta E1")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->beta_E1 = atof(tokptr);
			}
			else if(strcmp(tokptr, "timestep size")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->timestep_size = atof(tokptr);
			}
			else if(strcmp(tokptr, "input time")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->input_time = atof(tokptr);
			}
			else if(strcmp(tokptr, "external inhibition")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->extern_inhib = atof(tokptr);
			}
			else if(strcmp(tokptr, "global inhibition")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->global_inhibition = atof(tokptr);
			}
			else if(strcmp(tokptr, "vis1 strength")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->vis1_str = atof(tokptr);
			}
			else if(strcmp(tokptr, "vis1 location")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->vis1_loc = atof(tokptr);
			}
			else if(strcmp(tokptr, "vis1 sigma")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->vis1_sigma = atof(tokptr);
			}
			else if(strcmp(tokptr, "RC weight width")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->sigma_RC = atof(tokptr);
			}
			else if(strcmp(tokptr, "visring weight width")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->sigma_visring = atof(tokptr);
			}
			else if(strcmp(tokptr, "phi visring")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->phi_visring = atof(tokptr);
			}
			else if(strcmp(tokptr, "visring strength")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->visring_str = atof(tokptr);
			}			
			else if(strcmp(tokptr, "visring location")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->visring_loc = atof(tokptr);
			}
			else if(strcmp(tokptr, "visring sigma")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->visring_sigma = atof(tokptr);
			}
			else if(strcmp(tokptr, "learning rate")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->lrate = atof(tokptr);
			}
			else if(strcmp(tokptr, "learning")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->learning = atoi(tokptr);
			}
			else if(strcmp(tokptr, "normalise")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->normalise = atoi(tokptr);
			}
			else if(strcmp(tokptr, "LTD")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->LTD = atoi(tokptr);
			}
			else if(strcmp(tokptr, "dump length")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->weight_dump_length = atoi(tokptr);
			}
				
			else
				tokptr = strtok(NULL, "=");
			
		}
	}
	
	fclose(fptr);
	
	return;
}


void build_network(data *dptr, params *pptr, connect *cptr)
{
	int idx;
	
	pptr->timesteps = (int)(pptr->time/pptr->timestep_size);
	pptr->input_timesteps = (int)(pptr->input_time*(float)pptr->timesteps);
	
	dptr->time = (float *)malloc(pptr->timesteps * sizeof(float));
	
	dptr->activations_E1 = (float *)malloc(pptr->num_E1_cells*sizeof(float));
	dptr->prev_activations_E1 = (float *)malloc(pptr->num_E1_cells*sizeof(float));
	
	dptr->rates_E1 = (float *)malloc(pptr->num_E1_cells*sizeof(float));
	dptr->prev_rates_E1 = (float *)malloc(pptr->num_E1_cells*sizeof(float));
	
	dptr->rates_E1_time = (float **)malloc(pptr->num_E1_cells*sizeof(float *));
	dptr->rates_E1_time[0] = (float *)malloc((pptr->num_E1_cells*(pptr->timesteps/100))*sizeof(float));
	
	for(idx=1; idx<pptr->num_E1_cells; idx++)
	{
		dptr->rates_E1_time[idx] = dptr->rates_E1_time[0] + (idx * (pptr->timesteps/100));
	}
	
	dptr->activations_E1_time = (float **)malloc(pptr->num_E1_cells*sizeof(float *));
	dptr->activations_E1_time[0] = (float *)malloc((pptr->num_E1_cells*(pptr->timesteps/100))*sizeof(float));
	
	for(idx=1; idx<pptr->num_E1_cells; idx++)
	{
		dptr->activations_E1_time[idx] = dptr->activations_E1_time[0] + (idx * (pptr->timesteps/100));
	}
	
	cptr->RC1_connections = (int **)malloc(pptr->num_E1_cells*sizeof(int *));
	cptr->RC1_connections[0] = (int *)malloc((pptr->num_E1_cells*cptr->num_RC1_connections)*sizeof(int));
	
	for(idx=1; idx<pptr->num_E1_cells; idx++)
	{
		cptr->RC1_connections[idx] = cptr->RC1_connections[0] + (idx * cptr->num_RC1_connections);
	}
	
	cptr->RC1_weights =(float **)malloc(pptr->num_E1_cells*sizeof(float *));
	cptr->RC1_weights[0] = (float *)malloc((pptr->num_E1_cells*cptr->num_RC1_connections)*sizeof(float));
	
	for(idx=1; idx<pptr->num_E1_cells; idx++)
	{
		cptr->RC1_weights[idx] = cptr->RC1_weights[0] + (idx * cptr->num_RC1_connections);
	}
    
    dptr->favoured_view = (float *)malloc(pptr->num_E1_cells*sizeof(float));
    
    dptr->visual_input1 = (float *)malloc(pptr->num_E1_cells*sizeof(float));
	
	dptr->vis1_location_time = (float **)malloc(pptr->num_E1_cells*sizeof(float *));
	dptr->vis1_location_time[0] = (float *)malloc((pptr->num_E1_cells*(pptr->timesteps/100))*sizeof(float));
	
	for(idx=1; idx<pptr->num_E1_cells; idx++)
	{
		dptr->vis1_location_time[idx] = dptr->vis1_location_time[0] + (idx * (pptr->timesteps/100));
	}
	
	cptr->weights_visring =(float **)malloc(pptr->num_E1_cells*sizeof(float *));
	cptr->weights_visring[0] = (float *)malloc((pptr->num_E1_cells*pptr->num_E1_cells)*sizeof(float));
	
	for(idx=1; idx<pptr->num_E1_cells; idx++)
	{
		cptr->weights_visring[idx] = cptr->weights_visring[0] + (idx * pptr->num_E1_cells);
	}
	
	cptr->prev_weights_visring =(float **)malloc(pptr->num_E1_cells*sizeof(float *));
	cptr->prev_weights_visring[0] = (float *)malloc((pptr->num_E1_cells*pptr->num_E1_cells)*sizeof(float));
	
	for(idx=1; idx<pptr->num_E1_cells; idx++)
	{
		cptr->prev_weights_visring[idx] = cptr->prev_weights_visring[0] + (idx * pptr->num_E1_cells);
	}
	
	
	dptr->rates_visring = (float *)malloc(pptr->num_E1_cells*sizeof(float));
	dptr->prev_rates_visring = (float *)malloc(pptr->num_E1_cells*sizeof(float));

	dptr->visring_location_time = (float **)malloc(pptr->num_E1_cells*sizeof(float *));
	dptr->visring_location_time[0] = (float *)malloc((pptr->num_E1_cells*(pptr->timesteps/100))*sizeof(float));
	
	for(idx=1; idx<pptr->num_E1_cells; idx++)
	{
		dptr->visring_location_time[idx] = dptr->visring_location_time[0] + (idx * (pptr->timesteps/100));
	}
	
	dptr->E1_rates_average = (float *)malloc(pptr->num_E1_cells*sizeof(float));
	
	dptr->E1_rates_cumul = (float *)malloc(pptr->num_E1_cells*sizeof(float));
	
	dptr->i_total = (float*)malloc(pptr->num_E1_cells*sizeof(float));
		
	return;
}


void set_connectivity(connect *cptr, params *pptr)
{
	int cell, connection;
	
	for(cell=0; cell<pptr->num_E1_cells; cell++)
	{
		for(connection=0; connection<cptr->num_RC1_connections; connection++)
		{
			cptr->RC1_connections[cell][connection] = connection;
		}
	}
	
	return;
}


void set_favoured_view(data *dptr, params *pptr)
{
    int cell;
    float increment;
    
    increment = 360.0/(float)pptr->num_E1_cells;
    
    for(cell=0; cell<pptr->num_E1_cells; cell++)
    {
        dptr->favoured_view[cell] = (float)cell * increment;
    }
    
    return;
}


void set_RC_weights(params *pptr, data *dptr, connect *cptr)
{
	int cell, connection;
	float distance, distance1, distance2;
	float increment;
	FILE *fptr;
	
	increment = 360.0/(float)cptr->num_RC1_connections;
	
	for(cell=0; cell<pptr->num_E1_cells; cell++)
	{
		dptr->sumsq = 0.0;
		
		for(connection=0; connection<cptr->num_RC1_connections; connection++)
		{
			distance1 = fabs(((float)connection*increment) - dptr->favoured_view[cell]);
			distance2 = fabs(360.0-distance1);
			
			if(distance1<=distance2)
				distance = distance1;
			else
				distance = distance2;
			
			cptr->RC1_weights[cell][connection] = exp(-0.5 * (distance/pptr->sigma_RC) * (distance/pptr->sigma_RC));
			dptr->sumsq += cptr->RC1_weights[cell][connection] * cptr->RC1_weights[cell][connection];
			
		}
		dptr->sumsq = sqrt(dptr->sumsq);
		if(dptr->sumsq==0.0)
			continue;
		
		for(connection=0; connection<pptr->num_E1_cells; connection++)
		{
			cptr->RC1_weights[cell][connection] = cptr->RC1_weights[cell][connection] / dptr->sumsq;
		}
		
		
	}
	
	
	fptr = fopen("recurrent_weights.bdat", "wb");
	
	fwrite(&(cptr->RC1_weights[0][0]), sizeof(float), pptr->num_E1_cells*cptr->num_RC1_connections, fptr);
	
	fclose(fptr);
	
	return;
}

void set_visring_weights(params *pptr, data *dptr, connect *cptr)
{
	int cell, connection;
	float distance, distance1, distance2;
	float increment;
	
	increment = 360.0/(float)pptr->num_E1_cells;
	
	for(cell=0; cell<pptr->num_E1_cells; cell++)
	{
		dptr->sumsq = 0.0;
		
		for(connection=0; connection<pptr->num_E1_cells; connection++)
		{
			distance1 = fabs(((float)connection*increment) - dptr->favoured_view[cell]);
			distance2 = fabs(360.0-distance1);
			
			if(distance1<=distance2)
				distance = distance1;
			else
				distance = distance2;
			
			cptr->weights_visring[cell][connection] = exp(-0.5 * (distance/pptr->sigma_visring) * (distance/pptr->sigma_visring));
			dptr->sumsq += cptr->weights_visring[cell][connection] * cptr->weights_visring[cell][connection];
			
		}
		
		dptr->sumsq = sqrt(dptr->sumsq);
		if(dptr->sumsq==0.0)
			continue;
		
		for(connection=0; connection<pptr->num_E1_cells; connection++)
		{
			cptr->weights_visring[cell][connection] = cptr->weights_visring[cell][connection] / dptr->sumsq;
		}
		
		
	}
	return;
}



void free_memory(data *dptr, connect *cptr)
{
	free(dptr->activations_E1);
	free(dptr->prev_activations_E1);
	
	free(dptr->rates_E1);
	free(dptr->prev_rates_E1);
		
	free(dptr->rates_E1_time[0]);
	free(dptr->rates_E1_time);
	
	free(dptr->activations_E1_time[0]);
	free(dptr->activations_E1_time);
	    
    free(dptr->favoured_view);
    free(dptr->visual_input1);
	
	free(cptr->RC1_connections[0]);
	free(cptr->RC1_connections);
	
	free(cptr->RC1_weights[0]);
	free(cptr->RC1_weights);
	
	free(dptr->time);
	
	free(dptr->vis1_location_time[0]);
	free(dptr->vis1_location_time);
	
	free(cptr->weights_visring[0]);
	free(cptr->weights_visring);
	
	free(cptr->prev_weights_visring[0]);
	free(cptr->prev_weights_visring);
	
	free(dptr->rates_visring);
	free(dptr->prev_rates_visring);
	
	free(dptr->visring_location_time[0]);
	free(dptr->visring_location_time);
	
	free(dptr->E1_rates_average);
	free(dptr->E1_rates_cumul);
	
	free(dptr->i_total);

		
	return;
}
