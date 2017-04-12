#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_randist.h>

#include "init.h"
#include "data.h"
#include "file_io.h"
#include "array_utils.h"
#include "utils.h"
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
			
			if(strcmp(tokptr, "time")==0)			//Total time = [rotation time + non off time * num_HD_cells] * training_loops + pause time*2 + test time
			{
				tokptr = strtok(NULL, "=");
				pptr->time = atof(tokptr);
			}
			
			else if(strcmp(tokptr, "pause time")==0)			//Pause time during testing
			{
				tokptr = strtok(NULL, "=");
				pptr->pausetime = atof(tokptr);
			}
			else if(strcmp(tokptr, "rotation time")==0)			//Rotation time during training
			{
				tokptr = strtok(NULL, "=");
				pptr->rotation_time = atof(tokptr);
			}
			else if(strcmp(tokptr, "non offset training time")==0)		//This is the time the network waits at each location during training.
			{
				tokptr = strtok(NULL, "=");
				pptr->NonOffTrainTime = atof(tokptr);
			}
			else if(strcmp(tokptr, "training loops")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->training_loops = atoi(tokptr);			//Number of time to repeat total rotation time and total non-rotation time during training
			}			
			else if(strcmp(tokptr, "test time")==0)			//Rotation time during testing
			{
				tokptr = strtok(NULL, "=");
				pptr->test_time = atof(tokptr);
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
			else if(strcmp(tokptr, "HD to COMB connections")==0)
			{
				tokptr = strtok(NULL, "=");
				cptr->num_HD_COMB_connections = atoi(tokptr);
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
				pptr->numThreads = atoi(tokptr);
			}
			else if(strcmp(tokptr, "learning rate")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->lrate = atof(tokptr);
			}
			else if(strcmp(tokptr, "intro test")==0)
			{
				tokptr = strtok(NULL, "=");
				pptr->introtest = atoi(tokptr);
			}
			
			
			
				
			else
				tokptr = strtok(NULL, "=");
			
		}
	}
	
	fclose(fptr);
	
	return;
}


void build_network(data *dptr, params *pptr, connect *cptr, cArray **combptr, cArray **hdptr)
{
	//int idx;
	int connection;
	
		
	pptr->timesteps = (int)(pptr->time/pptr->timestep_size);
	pptr->test_timesteps = (int)(pptr->test_time/pptr->timestep_size);
	pptr->pausetimesteps = (int)(pptr->pausetime/pptr->timestep_size);
	pptr->rotation_timesteps = (int)(pptr->rotation_time/pptr->timestep_size);
	pptr->NonOffTrainTimesteps = (int)(pptr->NonOffTrainTime/pptr->timestep_size);
	
	
		
	pptr->conduction_buffer_size = (int)(pptr->conduction_delay/pptr->timestep_size);
	
	for (connection=0; connection<pptr->num_COMB_cells; connection++)
	{
		combptr[connection]->array = (float *)malloc(pptr->conduction_buffer_size *sizeof(float));
	}
	
	
	for (connection=0; connection<pptr->num_HD_cells; connection++)			
	{
		hdptr[connection]->array = (float *)malloc(pptr->conduction_buffer_size *sizeof(float)); 
	}
	
	dptr->time = (float *)malloc(pptr->timesteps * sizeof(float));
	
	dptr->activations_HD = (float *)malloc(pptr->num_HD_cells*sizeof(float));
	dptr->prev_activations_HD = (float *)malloc(pptr->num_HD_cells*sizeof(float));
	
	dptr->rates_HD = (float *)malloc(pptr->num_HD_cells*sizeof(float));
	dptr->prev_rates_HD = (float *)malloc(pptr->num_HD_cells*sizeof(float));
	
	
	dptr->rates_HD_time = get_2D_farray(pptr->num_HD_cells, (pptr->timesteps/100), 0);
	
	/*dptr->rates_HD_time = (float **)malloc(pptr->num_HD_cells*sizeof(float *));
	dptr->rates_HD_time[0] = (float *)malloc((pptr->num_HD_cells*(pptr->timesteps/100))*sizeof(float));
	
	for(idx=1; idx<pptr->num_HD_cells; idx++)
	{
		dptr->rates_HD_time[idx] = dptr->rates_HD_time[0] + (idx * (pptr->timesteps/100));
	}*/
	
	
	dptr->activations_HD_time = get_2D_farray(pptr->num_HD_cells, (pptr->timesteps/100), 0);
	
	/*dptr->activations_HD_time = (float **)malloc(pptr->num_HD_cells*sizeof(float *));
	dptr->activations_HD_time[0] = (float *)malloc((pptr->num_HD_cells*(pptr->timesteps/100))*sizeof(float));
	
	for(idx=1; idx<pptr->num_HD_cells; idx++)
	{
		dptr->activations_HD_time[idx] = dptr->activations_HD_time[0] + (idx * (pptr->timesteps/100));
	}*/	
	
	dptr->activations_COMB = (float *)malloc(pptr->num_COMB_cells*sizeof(float));
	dptr->prev_activations_COMB = (float *)malloc(pptr->num_COMB_cells*sizeof(float));
	
	dptr->rates_COMB = (float *)malloc(pptr->num_COMB_cells*sizeof(float));
	dptr->prev_rates_COMB = (float *)malloc(pptr->num_COMB_cells*sizeof(float));
	
	dptr->rates_COMB_time = get_2D_farray(pptr->num_COMB_cells, (pptr->timesteps/100), 0);

	/*dptr->rates_COMB_time = (float **)malloc(pptr->num_COMB_cells*sizeof(float *));
	dptr->rates_COMB_time[0] = (float *)malloc((pptr->num_COMB_cells*(pptr->timesteps/100))*sizeof(float));
	
	for(idx=1; idx<pptr->num_COMB_cells; idx++)
	{
		dptr->rates_COMB_time[idx] = dptr->rates_COMB_time[0] + (idx * (pptr->timesteps/100));
	}*/
	
	dptr->activations_COMB_time = get_2D_farray(pptr->num_COMB_cells, (pptr->timesteps/100), 0);

	/*dptr->activations_COMB_time = (float **)malloc(pptr->num_COMB_cells*sizeof(float *));
	dptr->activations_COMB_time[0] = (float *)malloc((pptr->num_COMB_cells*(pptr->timesteps/100))*sizeof(float));
	
	for(idx=1; idx<pptr->num_COMB_cells; idx++)
	{
		dptr->activations_COMB_time[idx] = dptr->activations_COMB_time[0] + (idx * (pptr->timesteps/100));
	}*/	
	
	
	dptr->rates_ROT = (float *)malloc(pptr->num_ROT_cells*sizeof(float));
	dptr->prev_rates_ROT = (float *)malloc(pptr->num_ROT_cells*sizeof(float));
	
	
	//dptr->rates_ROT_time = (float *)malloc((pptr->timesteps/100)*sizeof(float ));
		
	dptr->rates_ROT_time = get_2D_farray(pptr->num_ROT_cells, (pptr->timesteps/100), 0);
	
    dptr->favoured_view = (float *)malloc(pptr->num_HD_cells*sizeof(float));
    
    dptr->initialiser = (float *)malloc(pptr->num_HD_cells*sizeof(float));
	
	
	dptr->initialiser_location_time = get_2D_farray(pptr->num_HD_cells, (pptr->timesteps/100), 0);
	/*dptr->initialiser_location_time = (float **)malloc(pptr->num_HD_cells*sizeof(float *));
	dptr->initialiser_location_time[0] = (float *)malloc((pptr->num_HD_cells*(pptr->timesteps/100))*sizeof(float));
	
	for(idx=1; idx<pptr->num_HD_cells; idx++)
	{
		dptr->initialiser_location_time[idx] = dptr->initialiser_location_time[0] + (idx * (pptr->timesteps/100));
	}*/
		
    //Here's a matrix for the full possible HDtoCOMB connectivity, not the sample from sparse connectivity
	
	cptr->full_HDtoCOMB = get_2D_farray(pptr->num_COMB_cells,pptr->num_HD_cells, 0);
	/*cptr->full_HDtoCOMB =(float **)malloc(pptr->num_COMB_cells*sizeof(float *));
	cptr->full_HDtoCOMB[0] = (float *)malloc((pptr->num_COMB_cells*pptr->num_HD_cells)*sizeof(float));
	
	for(idx=1; idx<pptr->num_COMB_cells; idx++)
	{
		cptr->full_HDtoCOMB[idx] = cptr->full_HDtoCOMB[0] + (idx * pptr->num_HD_cells);
	}*/
	
		
	//cptr->ROT_COMB_weights =(float *)malloc(pptr->num_COMB_cells*sizeof(float ));
	//cptr->prev_ROT_COMB_weights =(float *)malloc(pptr->num_COMB_cells*sizeof(float ));
	
	cptr->ROT_COMB_weights = get_2D_farray(pptr->num_COMB_cells, pptr->num_ROT_cells, 0);
	cptr->prev_ROT_COMB_weights = get_2D_farray(pptr->num_COMB_cells, pptr->num_ROT_cells, 0);
			
			
	
	cptr->HD_COMB_connections = get_2D_iarray(pptr->num_COMB_cells, cptr->num_HD_COMB_connections, 0);
	
	 /*cptr->HD_COMB_connections = (int **)malloc(pptr->num_COMB_cells*sizeof(int *));
	 cptr->HD_COMB_connections[0] = (int *)malloc((cptr->num_HD_COMB_connections*pptr->num_COMB_cells)*sizeof(int));		//num_hd not num_hd_comb_connections??
	 
	 for(idx=1; idx<pptr->num_COMB_cells; idx++)
	 {
	 cptr->HD_COMB_connections[idx] = cptr->HD_COMB_connections[0] + (idx * cptr->num_HD_COMB_connections);
	 }*/
	 
	
	cptr->HD_COMB_weights = get_2D_farray(pptr->num_COMB_cells, cptr->num_HD_COMB_connections, 0);
	
	/*cptr->HD_COMB_weights =(float **)malloc(pptr->num_COMB_cells*sizeof(float *));
	cptr->HD_COMB_weights[0] = (float *)malloc((cptr->num_HD_COMB_connections*pptr->num_COMB_cells)*sizeof(float));
	
	for(idx=1; idx<pptr->num_COMB_cells; idx++)
	{
		cptr->HD_COMB_weights[idx] = cptr->HD_COMB_weights[0] + (idx * cptr->num_HD_COMB_connections);
	}*/
	
	
	cptr->prev_HD_COMB_weights = get_2D_farray(pptr->num_COMB_cells, cptr->num_HD_COMB_connections, 0);
	
	/*cptr->prev_HD_COMB_weights =(float **)malloc(pptr->num_COMB_cells*sizeof(float *));
	cptr->prev_HD_COMB_weights[0] = (float *)malloc((cptr->num_HD_COMB_connections*pptr->num_COMB_cells)*sizeof(float)); //Have I allocated this correctly?
	
	for(idx=1; idx<pptr->num_COMB_cells; idx++)
	{
		cptr->prev_HD_COMB_weights[idx] = cptr->prev_HD_COMB_weights[0] + (idx * cptr->num_HD_COMB_connections);
	}*/	

	
	cptr->COMB_HD_weights = get_2D_farray(pptr->num_HD_cells,pptr->num_COMB_cells,  0);
	
	/*cptr->COMB_HD_weights = (float **)malloc(pptr->num_HD_cells*sizeof(float *));
	cptr->COMB_HD_weights[0] = (float *)malloc((pptr->num_COMB_cells*pptr->num_HD_cells)*sizeof(float));
	
	for(idx=1; idx<pptr->num_HD_cells; idx++)
	{
		cptr->COMB_HD_weights[idx] = cptr->COMB_HD_weights[0] + (idx * pptr->num_COMB_cells);
	}*/		
	
	
	cptr->prev_COMB_HD_weights = get_2D_farray(pptr->num_HD_cells,pptr->num_COMB_cells, 0);
	
	
	/*cptr->prev_COMB_HD_weights =(float **)malloc(pptr->num_HD_cells*sizeof(float *));
	cptr->prev_COMB_HD_weights[0] = (float *)malloc((pptr->num_HD_cells*pptr->num_COMB_cells)*sizeof(float)); //Have I allocated this correctly?
	
	for(idx=1; idx<pptr->num_HD_cells; idx++)
	{
		cptr->prev_COMB_HD_weights[idx] = cptr->prev_COMB_HD_weights[0] + (idx * pptr->num_COMB_cells);
	}*/	

	

	return;
}


void set_connectivity(connect *cptr, params *pptr)
{
	int cell, connection;
	FILE *fptr;

	//HERE FOLLOWS SETUP OF RANDOM HD TO COMB CONNECTIVITY
	
	const gsl_rng_type *t=NULL;
	gsl_rng *r=NULL;

	gsl_rng_env_setup();
    
    t=gsl_rng_default;
    r=gsl_rng_alloc(t);
    
    gsl_rng_set(r, 1);
	
	int * a = malloc(pptr->num_HD_cells * sizeof(int));
	
	for(cell=0; cell<pptr->num_HD_cells; cell++)
	{
		a[cell] = cell;
	}
	
	
	
	
	for(cell=0; cell<pptr->num_COMB_cells; cell++)
	{
		gsl_ran_choose(r, &(cptr->HD_COMB_connections[cell][0]), cptr->num_HD_COMB_connections, a, pptr->num_HD_cells, sizeof(int));
	}					//(const gsl_rng * r, void * dest, size_t k, void * src, size_t n, size_t size)
		
	gsl_rng_free(r);
	
	free(a);

	 
	 
	fptr = file_open("HD_COMB_connectivity.txt", "w");
	for(cell=0;cell<pptr->num_COMB_cells;cell++)
	{
		for(connection=0;connection<cptr->num_HD_COMB_connections;connection++)
		{
				fprintf(fptr, "%d ", cptr->HD_COMB_connections[cell][connection]);
		}
		fprintf(fptr,"\n");
	}
		
	
	fclose(fptr);
	
	
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
	int cell, connection;
	
	srand48(time(NULL));


		
	for(cell=0; cell<pptr->num_COMB_cells; cell++)
	{
		for(connection=0; connection<cptr->num_HD_COMB_connections; connection++)
		{
			cptr->HD_COMB_weights[cell][connection] = drand48();//((float)rand() / (float)RAND_MAX)* 0.00001;
		}
	}

	normalise_HDtoCOMB_weights(cptr, pptr);
	return;
	
	
}

void set_ROT_COMB_weights(params *pptr, data *dptr, connect *cptr)
{
	int cell;
	int connection;
	srand48(time(NULL));

	
	for(cell=0; cell<pptr->num_COMB_cells;cell++)
	{
		for(connection=0; connection<pptr->num_ROT_cells; connection++)
		{
			cptr->ROT_COMB_weights[cell][connection] = drand48();//((float)rand() / (float)RAND_MAX))* 0.1;
		}
	}
	
	normalise_ROTtoCOMB_weights(cptr, pptr);
	
	return;
}


void set_COMB_HD_weights(params *pptr, data *dptr, connect *cptr)
{

	int cell, connection;
	srand48(time(NULL));

	
	for(cell=0; cell<pptr->num_HD_cells; cell++)
	{
		for(connection=0; connection<pptr->num_COMB_cells; connection++)
		{
			cptr->COMB_HD_weights[cell][connection] = drand48(); //(float)rand() / (float)RAND_MAX)* 0.1;
		}
	}
	
	normalise_COMBtoHD_weights(cptr, pptr);
	
	
	return;
	
	
	
}


void free_memory(params *pptr, data *dptr, connect *cptr, cArray **combptr, cArray **hdptr)
{
	int connection;
	
	free(dptr->time);
		

	for(connection=0; connection<pptr->num_COMB_cells; connection++)
	{
		free(combptr[connection]->array);
	}
	
	for(connection=0; connection<pptr->num_HD_cells; connection++)
	{
		free(hdptr[connection]->array);
	}
	
	free(dptr->activations_HD);
	free(dptr->prev_activations_HD);
	
	free(dptr->rates_HD);
	free(dptr->prev_rates_HD);
	
	free_2D_farray(dptr->rates_HD_time);
	/*free(dptr->rates_HD_time[0]);
	free(dptr->rates_HD_time);*/
	
	free_2D_farray(dptr->activations_HD_time);
	/*free(dptr->activations_HD_time[0]);
	free(dptr->activations_HD_time);*/
	
	free(dptr->activations_COMB);
	free(dptr->prev_activations_COMB);
	
	free(dptr->rates_COMB);
	free(dptr->prev_rates_COMB);
	
	free_2D_farray(dptr->rates_COMB_time);
	free_2D_farray(dptr->activations_COMB_time);
	/*free(dptr->rates_COMB_time[0]);
	free(dptr->rates_COMB_time);
	
	free(dptr->activations_COMB_time[0]);
	free(dptr->activations_COMB_time);*/
	
	free(dptr->rates_ROT);
	free(dptr->prev_rates_ROT);
	
	free_2D_farray(dptr->rates_ROT_time);
	
	free_2D_iarray(cptr->HD_COMB_connections);
	/*free(cptr->HD_COMB_connections[0]);
	free(cptr->HD_COMB_connections);*/

	
	/*free(cptr->HD_COMB_weights[0]);
	free(cptr->HD_COMB_weights);
		
	free(cptr->prev_HD_COMB_weights[0]);
	free(cptr->prev_HD_COMB_weights);
		
	free(cptr->COMB_HD_weights[0]);
	free(cptr->COMB_HD_weights);
	
	free(cptr->prev_COMB_HD_weights[0]);
	free(cptr->prev_COMB_HD_weights);*/
	
	free_2D_farray(cptr->HD_COMB_weights);
	
	free_2D_farray(cptr->COMB_HD_weights);
	
	free_2D_farray(cptr->prev_HD_COMB_weights); 	
	free_2D_farray(cptr->prev_COMB_HD_weights);

	
	

	free(dptr->favoured_view);
    free(dptr->initialiser);
	
	free_2D_farray(dptr->initialiser_location_time);
	/*free(dptr->initialiser_location_time[0]);
	free(dptr->initialiser_location_time);*/
	
	free_2D_farray(cptr->full_HDtoCOMB);
	/*free(cptr->full_HDtoCOMB[0]);
	free(cptr->full_HDtoCOMB);*/

	/*free(cptr->ROT_COMB_weights[0]);
	free(cptr->ROT_COMB_weights);
	
	free(cptr->prev_ROT_COMB_weights[0]);
	free(cptr->prev_ROT_COMB_weights);
	*/
	
	free_2D_farray(cptr->ROT_COMB_weights);
	free_2D_farray(cptr->prev_ROT_COMB_weights);
		

	return;
}
