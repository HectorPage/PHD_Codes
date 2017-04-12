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
			else if(strcmp(tokptr, "HD cells")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->num_HD_cells = atoi(tokptr);
			}
			else if(strcmp(tokptr, "COMB cells")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->num_COMB_cells = atoi(tokptr);
			}
			else if(strcmp(tokptr, "ROT cells")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->num_ROT_cells = atoi(tokptr);
			}
			else if(strcmp(tokptr, "NOROT cells")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->num_NOROT_cells = atoi(tokptr);
			}
			else if(strcmp(tokptr, "phi HDCOMB")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->phi_HDCOMB = atof(tokptr);
			}
			else if(strcmp(tokptr, "phi COMBHD")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->phi_COMBHD = atof(tokptr);
			}
			else if(strcmp(tokptr, "phi ROT")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->phi_ROT = atof(tokptr);
			}
			else if(strcmp(tokptr, "phi NOROT")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->phi_NOROT = atof(tokptr);
			}			
			else if(strcmp(tokptr, "tau HD")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->tau_HD = atof(tokptr);
			}
			else if(strcmp(tokptr, "tau COMB")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->tau_COMB = atof(tokptr);
			}
			else if(strcmp(tokptr, "alpha HD")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->alpha_HD = atof(tokptr);
			}
			else if(strcmp(tokptr, "beta HD")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->beta_HD = atof(tokptr);
			}
			else if(strcmp(tokptr, "alpha COMB")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->alpha_COMB = atof(tokptr);
			}
			else if(strcmp(tokptr, "beta COMB")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->beta_COMB = atof(tokptr);
			}
			else if(strcmp(tokptr, "timestep size")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->timestep_size = atof(tokptr);
			}
			else if(strcmp(tokptr, "pause time")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->pausetime = atof(tokptr);
			}
			else if(strcmp(tokptr, "rotation time")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->rotation_time = atof(tokptr);
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
			else if(strcmp(tokptr, "COMB inhibition")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->COMB_inhibition = atof(tokptr);
			}
			else if(strcmp(tokptr, "initialiser strength")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->initialiser_str = atof(tokptr);
			}
			else if(strcmp(tokptr, "initialiser location")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->initialiser_loc = atof(tokptr);
			}
			else if(strcmp(tokptr, "initialiser sigma")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->initialiser_sigma = atof(tokptr);
			}
			else if(strcmp(tokptr, "HD to COMB weight width")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->sigma_HD_COMB = atof(tokptr);
			}
			else if(strcmp(tokptr, "COMB to HD weight width")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->sigma_COMB_HD = atof(tokptr);
			}
			else if(strcmp(tokptr, "conduction delay")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->conduction_delay = atof(tokptr);
			}			
			else if(strcmp(tokptr, "velocity")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->velocity = atof(tokptr);
			}
			else if(strcmp(tokptr, "threads")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->numThreads = atof(tokptr);
			}
			
				
			else
				tokptr = strtok(NULL, "=");
			
		}
	}
	
	fclose(fptr);
	
	return;
}


void build_network(data *dptr, params *pptr, connect *cptr, cArray **bptr, cArray **bxptr)
{
	int idx;
	int connection;
	
	pptr->timesteps = (int)(pptr->time/pptr->timestep_size);
	pptr->input_timesteps = (int)(pptr->input_time/pptr->timestep_size);
	pptr->pausetimesteps = (int)(pptr->pausetime/pptr->timestep_size);
	pptr->rotation_timesteps = (int)(pptr->rotation_time/pptr->timestep_size);
	
	
	pptr->conduction_buffer_size = (int)(pptr->conduction_delay/pptr->timestep_size);
	
	for (connection=0; connection<pptr->num_COMB_cells; connection++)
	{
		bptr[connection]->array = (float *)malloc(pptr->conduction_buffer_size *sizeof(float));
	}
	
	
	for (connection=0; connection<pptr->num_HD_cells; connection++)
	{
		bxptr[connection]->array = (float *)malloc(pptr->conduction_buffer_size *sizeof(float));
	}
	
	
	dptr->time = (float *)malloc(pptr->timesteps * sizeof(float));
	
	dptr->activations_HD = (float *)malloc(pptr->num_HD_cells*sizeof(float));
	dptr->prev_activations_HD = (float *)malloc(pptr->num_HD_cells*sizeof(float));
	
	dptr->rates_HD = (float *)malloc(pptr->num_HD_cells*sizeof(float));
	dptr->prev_rates_HD = (float *)malloc(pptr->num_HD_cells*sizeof(float));
	
	dptr->rates_HD_time = (float **)malloc(pptr->num_HD_cells*sizeof(float *));
	dptr->rates_HD_time[0] = (float *)malloc((pptr->num_HD_cells*(pptr->timesteps/100))*sizeof(float));
	
	for(idx=1; idx<pptr->num_HD_cells; idx++)
	{
		dptr->rates_HD_time[idx] = dptr->rates_HD_time[0] + (idx * (pptr->timesteps/100));
	}
	
	dptr->activations_HD_time = (float **)malloc(pptr->num_HD_cells*sizeof(float *));
	dptr->activations_HD_time[0] = (float *)malloc((pptr->num_HD_cells*(pptr->timesteps/100))*sizeof(float));
	
	for(idx=1; idx<pptr->num_HD_cells; idx++)
	{
		dptr->activations_HD_time[idx] = dptr->activations_HD_time[0] + (idx * (pptr->timesteps/100));
	}	
	
	dptr->activations_COMB = (float *)malloc(pptr->num_COMB_cells*sizeof(float));
	dptr->prev_activations_COMB = (float *)malloc(pptr->num_COMB_cells*sizeof(float));
	
	dptr->rates_COMB = (float *)malloc(pptr->num_COMB_cells*sizeof(float));
	dptr->prev_rates_COMB = (float *)malloc(pptr->num_COMB_cells*sizeof(float));
	
	dptr->rates_COMB_time = (float **)malloc(pptr->num_COMB_cells*sizeof(float *));
	dptr->rates_COMB_time[0] = (float *)malloc((pptr->num_COMB_cells*(pptr->timesteps/100))*sizeof(float));
	
	for(idx=1; idx<pptr->num_COMB_cells; idx++)
	{
		dptr->rates_COMB_time[idx] = dptr->rates_COMB_time[0] + (idx * (pptr->timesteps/100));
	}
	
	dptr->activations_COMB_time = (float **)malloc(pptr->num_COMB_cells*sizeof(float *));
	dptr->activations_COMB_time[0] = (float *)malloc((pptr->num_COMB_cells*(pptr->timesteps/100))*sizeof(float));
	
	for(idx=1; idx<pptr->num_COMB_cells; idx++)
	{
		dptr->activations_COMB_time[idx] = dptr->activations_COMB_time[0] + (idx * (pptr->timesteps/100));
	}	
	
	dptr->rates_ROT_time = (float *)malloc((pptr->timesteps/100)*sizeof(float *));
	dptr->rates_NOROT_time = (float *)malloc((pptr->timesteps/100)*sizeof(float *));

	cptr->HD_COMB_connections = (int **)malloc(pptr->num_COMB_cells*sizeof(int *));
	cptr->HD_COMB_connections[0] = (int *)malloc((pptr->num_HD_cells*pptr->num_COMB_cells)*sizeof(int));
	
	for(idx=1; idx<pptr->num_COMB_cells; idx++)
	{
		cptr->HD_COMB_connections[idx] = cptr->HD_COMB_connections[0] + (idx * pptr->num_HD_cells);
	}
	
	cptr->HD_COMB_weights =(float **)malloc(pptr->num_COMB_cells*sizeof(float *));
	cptr->HD_COMB_weights[0] = (float *)malloc((pptr->num_HD_cells*pptr->num_COMB_cells)*sizeof(float));
	
	for(idx=1; idx<pptr->num_COMB_cells; idx++)
	{
		cptr->HD_COMB_weights[idx] = cptr->HD_COMB_weights[0] + (idx * pptr->num_HD_cells);
	}
    
	cptr->COMB_HD_connections = (int **)malloc(pptr->num_HD_cells*sizeof(int *));
	cptr->COMB_HD_connections[0] = (int *)malloc((pptr->num_COMB_cells*pptr->num_HD_cells)*sizeof(int)); //swapped the nums around to help with seg faults?
	
	for(idx=1; idx<pptr->num_HD_cells; idx++)
	{
		cptr->COMB_HD_connections[idx] = cptr->COMB_HD_connections[0] + (idx * pptr->num_COMB_cells);
	}
	
	cptr->COMB_HD_weights = (float **)malloc(pptr->num_HD_cells*sizeof(float *));
	cptr->COMB_HD_weights[0] = (float *)malloc((pptr->num_HD_cells*pptr->num_COMB_cells)*sizeof(float));
	
	for(idx=1; idx<pptr->num_HD_cells; idx++)
	{
		cptr->COMB_HD_weights[idx] = cptr->COMB_HD_weights[0] + (idx * pptr->num_COMB_cells);
	}	
	
    dptr->favoured_view = (float *)malloc(pptr->num_HD_cells*sizeof(float));
    
    dptr->initialiser = (float *)malloc(pptr->num_HD_cells*sizeof(float));
	
	dptr->initialiser_location_time = (float **)malloc(pptr->num_HD_cells*sizeof(float *));
	dptr->initialiser_location_time[0] = (float *)malloc((pptr->num_HD_cells*(pptr->timesteps))*sizeof(float));
	
	for(idx=1; idx<pptr->num_HD_cells; idx++)
	{
		dptr->initialiser_location_time[idx] = dptr->initialiser_location_time[0] + (idx * (pptr->timesteps));
	}

		
	return;
}


void set_connectivity(connect *cptr, params *pptr)
{
	int cell, connection;
	
	for(cell=0; cell<pptr->num_HD_cells; cell++)
	{
		for(connection=0; connection<pptr->num_COMB_cells; connection++)
		{
			cptr->COMB_HD_connections[cell][connection] = connection;
		}
	}
	
	for(cell=0; cell<pptr->num_COMB_cells; cell++)
	{
		for(connection=0; connection<pptr->num_HD_cells; connection++)
		{
			cptr->HD_COMB_connections[cell][connection] = connection;
		
		}
	}
		
	
	return;
}


void set_favoured_view(data *dptr, params *pptr)
{
    int cell;
    float increment;
    
    increment = 360.0/(float)pptr->num_HD_cells;
    
    for(cell=0; cell<pptr->num_HD_cells; cell++)
    {
        dptr->favoured_view[cell] = (float)cell * increment;
    }
    
    return;
}


void set_HD_COMB_weights(params *pptr, data *dptr, connect *cptr)
{
	int cell, connection, COMBfav;
	float distance, distance1, distance2;
	float new_direction;
	float offset;
	float increment;
	FILE *fptr;
	
	increment = 360.0/(float)pptr->num_HD_cells;	
	offset = (pptr->velocity * pptr->conduction_delay);		

	
	COMBfav = 0; //variable used to make sure COMB cells have correct favoured direction for each of the two COMB cell populations, rather than trying to get favouried_view for cell>#HD cells (which doesn't exist!)
	
	for(cell=0; cell<pptr->num_COMB_cells/2; cell++)
	{
		for(connection=0; connection<pptr->num_HD_cells; connection++)
		{
			
			new_direction = offset + dptr->favoured_view[connection];
			if(new_direction>360.0)
				new_direction = new_direction -360.0;
			
			distance1 = fabs(new_direction - dptr->favoured_view[COMBfav]);
			distance2 = fabs(360.0-distance1);
			
			if(distance1<=distance2)
				distance = distance1;
			else
				distance = distance2;
			
			cptr->HD_COMB_weights[cell][connection] = exp(-0.5 * (distance/pptr->sigma_HD_COMB) * (distance/pptr->sigma_HD_COMB));
			
			
		}
		
				
		COMBfav++;
	}
	
	COMBfav = 0;
	
	for(cell=pptr->num_COMB_cells/2; cell<pptr->num_COMB_cells; cell++)
	{		
		for(connection=0; connection<pptr->num_HD_cells; connection++)
		{
			distance1 = fabs(((float)connection*increment) - dptr->favoured_view[COMBfav]);
			distance2 = fabs(360.0-distance1);
			
			if(distance1<=distance2)
				distance = distance1;
			else
				distance = distance2;
			
			cptr->HD_COMB_weights[cell][connection] = exp(-0.5 * (distance/pptr->sigma_HD_COMB) * (distance/pptr->sigma_HD_COMB));
						
		}
		
	
		COMBfav++;
	}
	
	for(cell=0;cell<pptr->num_COMB_cells;cell++)
	{
		dptr->sumsq_HD_COMB = 0.0;	
		
		for(connection=0; connection<pptr->num_HD_cells; connection++)
		{
			dptr->sumsq_HD_COMB += cptr->HD_COMB_weights[cell][connection] * cptr->HD_COMB_weights[cell][connection];
			
		}
		
		dptr->sumsq_HD_COMB = sqrt(dptr->sumsq_HD_COMB);
		if(dptr->sumsq_HD_COMB==0.0)
			continue;
		
		for(connection=0; connection<pptr->num_HD_cells; connection++)
		{
			cptr->HD_COMB_weights[cell][connection] = cptr->HD_COMB_weights[cell][connection] / dptr->sumsq_HD_COMB;
		}
		
		
	}
	
	
	
	
	
	fptr = fopen("HDtoCOMBweights.bdat", "wb");
	
	fwrite(&(cptr->HD_COMB_weights[0][0]), sizeof(float), pptr->num_HD_cells*pptr->num_COMB_cells, fptr);
	
	fclose(fptr);
	
	return;
	
	
}


void set_COMB_HD_weights(params *pptr, data *dptr, connect *cptr)
{

	int cell, connection, COMBfav;
	float distance, distance1, distance2;
	float new_direction;
	float offset;
	FILE *fptr;
	
	offset = (pptr->velocity * pptr->conduction_delay);		
		
	for(cell=0; cell<pptr->num_HD_cells; cell++)
	{
		COMBfav = 0;	//variable used to make sure COMB cells have correct favoured direction for each of the two COMB cell populations, rather than trying to get favouried_view for cell>#HD cells (which doesn't exist!)
		
			for(connection=0; connection<pptr->num_COMB_cells/2; connection++)
			{

				new_direction = offset + dptr->favoured_view[COMBfav];
				if(new_direction>360.0)
				new_direction = new_direction -360.0;
			
				distance1 = fabs(dptr->favoured_view[cell]-new_direction);
				distance2 = fabs(360.0-distance1);
			
				if(distance1<=distance2)
					distance = distance1;
				else {
					distance = distance2;
				}
				cptr->COMB_HD_weights[cell][connection] = exp(-0.5 * (distance/pptr->sigma_COMB_HD) * (distance/pptr->sigma_COMB_HD));
				

				COMBfav++;
			}
	}

	
	for(cell=0; cell<pptr->num_HD_cells; cell++)
	{
		
		COMBfav = 0;

					
			for(connection=pptr->num_COMB_cells/2; connection<pptr->num_COMB_cells; connection++)
			{
				
				distance1 = fabs(dptr->favoured_view[cell] - dptr->favoured_view[COMBfav]);
				distance2 = fabs(360.0-distance1);
			
				if(distance1<=distance2)
					distance = distance1;
				else
					distance = distance2;
				
				cptr->COMB_HD_weights[cell][connection] = exp(-0.5 * (distance/pptr->sigma_COMB_HD) * (distance/pptr->sigma_COMB_HD));
								
				COMBfav++;
			}
		

			
	}
	
			
		
	for(cell=0;cell<pptr->num_HD_cells;cell++)
	{
		dptr->sumsq_COMB_HD = 0.0;	
						
		for(connection=0; connection<pptr->num_COMB_cells; connection++)
		{
			dptr->sumsq_COMB_HD += cptr->COMB_HD_weights[cell][connection] * cptr->COMB_HD_weights[cell][connection];
		
		}
		
		dptr->sumsq_COMB_HD = sqrt(dptr->sumsq_COMB_HD);
		if(dptr->sumsq_COMB_HD==0.0)
			continue;
		
		for(connection=0; connection<pptr->num_COMB_cells; connection++)
		{
			cptr->COMB_HD_weights[cell][connection] = cptr->COMB_HD_weights[cell][connection] / dptr->sumsq_COMB_HD;
		}


	}
	
	
	fptr = fopen("COMBtoHDweights.bdat", "wb");
	
	fwrite(&(cptr->COMB_HD_weights[0][0]), sizeof(float), pptr->num_HD_cells*pptr->num_COMB_cells, fptr);
	
	fclose(fptr);
	
	return;
	
	
}


void free_memory(params *pptr, data *dptr, connect *cptr, cArray **bptr, cArray **bxptr)
{
	int connection;
	
	free(dptr->activations_HD);
	free(dptr->prev_activations_HD);
	
	free(dptr->rates_HD);
	free(dptr->prev_rates_HD);
	
	free(dptr->activations_COMB);
	free(dptr->prev_activations_COMB);
	
	free(dptr->rates_COMB);
	free(dptr->prev_rates_COMB);
	
	free(dptr->rates_HD_time[0]);
	free(dptr->rates_HD_time);
	
	free(dptr->activations_HD_time[0]);
	free(dptr->activations_HD_time);
	
	free(dptr->rates_COMB_time[0]);
	free(dptr->rates_COMB_time);
	
	free(dptr->activations_COMB_time[0]);
	free(dptr->activations_COMB_time);
	
	free(dptr->rates_ROT_time);
	free(dptr->rates_NOROT_time);
	
	free(cptr->COMB_HD_weights[0]);
	free(cptr->COMB_HD_weights);
	free(cptr->COMB_HD_connections[0]);
	free(cptr->COMB_HD_connections);
	
	free(cptr->HD_COMB_weights[0]);
	free(cptr->HD_COMB_weights);
	free(cptr->HD_COMB_connections[0]);
	free(cptr->HD_COMB_connections);
    
    free(dptr->favoured_view);
    free(dptr->initialiser);
		
	free(dptr->time);
	
	free(dptr->initialiser_location_time[0]);
	free(dptr->initialiser_location_time);
	
	for(connection=0; connection<pptr->num_COMB_cells; connection++)
	{
		free(bptr[connection]->array);
	}
	
	for(connection=0; connection<pptr->num_HD_cells; connection++)
	{
		free(bxptr[connection]->array);
	}
	
	
		
	return;
}
