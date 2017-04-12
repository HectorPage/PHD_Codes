#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_randist.h>

#include "init.h"
#include "data.h"
#include "file_io.h"
#include "calculate.h"

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
			else if(strcmp(tokptr, "phi RC1 test")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->phi_RC1_test = atof(tokptr);
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
			else if(strcmp(tokptr, "training time")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->training_time = atof(tokptr);
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
			else if(strcmp(tokptr, "input strength")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->input_str = atof(tokptr);
			}
			else if(strcmp(tokptr, "input location")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->input_loc = atof(tokptr);
			}
			else if(strcmp(tokptr, "input sigma")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->input_sigma = atof(tokptr);
			}
			else if(strcmp(tokptr, "normalise")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->normalise = atoi(tokptr);
			}
			else if(strcmp(tokptr, "velocity")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->velocity = atof(tokptr);
			}
			else if(strcmp(tokptr, "conduction delay")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->conduction_delay = atof(tokptr);
			}			
			else if(strcmp(tokptr, "parallel threads")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->threads = atoi(tokptr);
			}
			else if(strcmp(tokptr, "sigmoid")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->sigmoid = atoi(tokptr);
			}
			else if(strcmp(tokptr, "symmetrical testing")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->symtest = atoi(tokptr);
			}
			else if(strcmp(tokptr, "symmetrical sigma")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->sym_sigma = atof(tokptr);
			}
			else if(strcmp(tokptr, "symmetrical strength")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->sym_str = atof(tokptr);
			}
			else
				tokptr = strtok(NULL, "=");
			
		}
	}
	
	fclose(fptr);
	
	return;
}


void build_network(data *dptr, params *pptr, connect *cptr, cArray **bptr)
{
	int idx;
	int connection;
	
	pptr->timesteps = (int)(pptr->time/pptr->timestep_size);
	pptr->training_timesteps = (int)(pptr->training_time/(float)pptr->timestep_size);
	
	//Here's the conduction buffer allocation stuff
	
	pptr->conduction_buffer_size = (int)(pptr->conduction_delay/pptr->timestep_size);
	
	for (connection=0; connection<pptr->num_E1_cells; connection++)						
	{
		bptr[connection]->array = (float *)malloc(pptr->conduction_buffer_size *sizeof(float));
	}
	
	
	
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
	
	//Here's a matrix for the full possible RC connectivity, not the sample from sparse connectivity
	
	cptr->full_RC =(float **)malloc(pptr->num_E1_cells*sizeof(float *));
	cptr->full_RC[0] = (float *)malloc((pptr->num_E1_cells*pptr->num_E1_cells)*sizeof(float));
	
	for(idx=1; idx<pptr->num_E1_cells; idx++)
	{
		cptr->full_RC[idx] = cptr->full_RC[0] + (idx * pptr->num_E1_cells);
	}
	

	cptr->prev_weights_RC1 =(float **)malloc(pptr->num_E1_cells*sizeof(float *));
	cptr->prev_weights_RC1[0] = (float *)malloc((pptr->num_E1_cells*cptr->num_RC1_connections)*sizeof(float));
	
	for(idx=1; idx<pptr->num_E1_cells; idx++)
	{
		cptr->prev_weights_RC1[idx] = cptr->prev_weights_RC1[0] + (idx * cptr->num_RC1_connections);
	}	
    
    dptr->favoured_view = (float *)malloc(pptr->num_E1_cells*sizeof(float));
    
    dptr->input = (float *)malloc(pptr->num_E1_cells*sizeof(float));
	
	dptr->input_location_time = (float **)malloc(pptr->num_E1_cells*sizeof(float *));
	dptr->input_location_time[0] = (float *)malloc((pptr->num_E1_cells*(pptr->timesteps/100))*sizeof(float));
	
	for(idx=1; idx<pptr->num_E1_cells; idx++)
	{
		dptr->input_location_time[idx] = dptr->input_location_time[0] + (idx * (pptr->timesteps/100));
	}
	
	//Here's stuff for LTD
	
	dptr->E1_rates_average = (float *)malloc(pptr->num_E1_cells*sizeof(float));
	
	dptr->E1_rates_cumul = (float *)malloc(pptr->num_E1_cells*sizeof(float));
	
	
	//Here's the stuff for recording and storing weight vectors at various points in simulation
	
	dptr->weight_vector_intermediate = (float *)malloc(pptr->num_E1_cells*sizeof(float));
	
	dptr->weight_vector_final = (float *)malloc(pptr->num_E1_cells*sizeof(float));
	
	dptr->weight_vector = (float *)malloc(pptr->num_E1_cells*sizeof(float));
	

	return;
}


void set_connectivity(connect *cptr, params *pptr)
{
	int cell;
	const gsl_rng_type *t=NULL;
    gsl_rng *r=NULL;
	
	gsl_rng_env_setup();
    
    t=gsl_rng_default;
    r=gsl_rng_alloc(t);
    
    gsl_rng_set(r, 1);
	
	int a[pptr->num_E1_cells];
	
	for(cell=0; cell<pptr->num_E1_cells; cell++)
	{
		a[cell] = cell;
	}
	
	
	
	for(cell=0; cell<pptr->num_E1_cells; cell++)
	{
		gsl_ran_choose(r, &(cptr->RC1_connections[cell][0]), cptr->num_RC1_connections, a, pptr->num_E1_cells, sizeof(int)); 

	}
	
	gsl_rng_free(r);
	
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
	int cell, connection, presynaptic;
	float increment, distance1, distance2, distance;
	float new_direction, offset;
	
	increment = 360.0/(float)cptr->num_RC1_connections;
	
	offset = (pptr->velocity * pptr->conduction_delay);
	
#pragma omp for private(cell, connection, distance1, distance2, distance, presynaptic)
	for(cell=0; cell<pptr->num_E1_cells; cell++)
	{
		
		for(connection=0; connection<cptr->num_RC1_connections; connection++)
		{
			
			presynaptic = cptr->RC1_connections[cell][connection];
			
			new_direction = dptr->favoured_view[presynaptic] + offset;    
			
			if(new_direction>360.0)
				new_direction = new_direction - 360.0;
			if(new_direction<=0)
				new_direction = new_direction + 360.0;
			
			
			
			distance1 = fabs(dptr->favoured_view[cell] - new_direction);			
			distance2 = fabs(360.0 - distance1);
			
			if(distance1<=distance2)
				distance = distance1;
			else
			{
				distance = distance2;
			}
			
			cptr->RC1_weights[cell][connection] = exp(-0.5 * (distance/pptr->sym_sigma) * (distance/pptr->sym_sigma));
			
		}
	}
	
		
	normalise_RC_weights(cptr, pptr);
	
	return;
}




void free_memory(data *dptr, connect *cptr, cArray **bptr, params *pptr)
{
	int connection;
	
	free(dptr->activations_E1);
	free(dptr->prev_activations_E1);

	free(dptr->rates_E1);
	free(dptr->prev_rates_E1);
			
	free(dptr->rates_E1_time[0]);
	free(dptr->rates_E1_time);
	
	free(dptr->activations_E1_time[0]);
	free(dptr->activations_E1_time);
	
    free(dptr->favoured_view);
    free(dptr->input);
		
	free(cptr->RC1_connections[0]);
	free(cptr->RC1_connections);
	
	free(cptr->RC1_weights[0]);
	free(cptr->RC1_weights);
		
	free(cptr->full_RC[0]);
	free(cptr->full_RC);
		
	free(cptr->prev_weights_RC1[0]);
	free(cptr->prev_weights_RC1);
		
	free(dptr->time);
		
	free(dptr->input_location_time[0]);
	free(dptr->input_location_time);
			
	for(connection=0; connection<pptr->num_E1_cells; connection++)
	{
		free(bptr[connection]->array);
	}
	
	free(dptr->E1_rates_average);
	free(dptr->E1_rates_cumul);
	
	free(dptr->weight_vector);
	free(dptr->weight_vector_intermediate);
	free(dptr->weight_vector_final);
		
	return;
}
