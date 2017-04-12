#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "calculate.h"
#include "data.h"
#include "file_io.h"
#include "init.h"



int main(int argc, char **argv)
{
	char *fname;
	int timestep;
	params *pptr;
	data *dptr;
	connect *cptr;
	cArray **bptr;
	int counter;
	float increment;
	int connection;
	double start, end;
	
	if(argc!=3)
	{
		fprintf(stderr, "Usage: toy -fname <parameters_file>\n");
		exit(1);
	}
	else
	{
		pptr = (params *)malloc(sizeof(params));
		cptr = (connect *)malloc(sizeof(connect));
		fname = argv[2];
	}
	
	read_parameters(fname, pptr, cptr);
	
	dptr = (data *)malloc(sizeof(data));
	
	bptr = (cArray **)malloc(sizeof(cArray *) * pptr->num_E1_cells);
	
	for(counter=0; counter<pptr->num_E1_cells; counter++) 
	{
		bptr[counter] = malloc(sizeof(cArray));
	}
	
	
	build_network(dptr, pptr, cptr, bptr);
	
	
	set_connectivity(cptr, pptr);
	
    set_favoured_view(dptr, pptr);
	
	set_RC_weights(pptr, dptr, cptr); //N.B. - These are intialised at zero, to speed up learning process.
	
		
	load_full_RC(cptr, pptr);
	
	
	save_RCweights_init(cptr, pptr);
	
	
	
	/* Simulate Network */
	
	zero_states(dptr, pptr, bptr, cptr);
	zero_input(dptr, pptr);

		
	printf("\nNetwork now being simulated...\n\n");
	fflush(stdout);
	
	
	printf("Training Phase...\n");
	fflush(stdout);

	
	omp_set_num_threads(pptr->threads);
	start = omp_get_wtime();
	
#pragma omp parallel private (timestep)
{
#pragma omp single
	{
		pptr->testphase = 0; 
	}
	
	for(timestep=0; timestep<pptr->start_length; timestep++)
	{
#pragma omp single nowait
		{
			printf("\rTimestep %d", timestep+1);
			fflush(stdout);
		}
		
		
		record_previous_states(dptr, pptr, cptr);
		
		
#pragma omp single
		{
			increment = pptr->velocity*pptr->timestep_size;
			pptr->input_loc += increment;
			if(pptr->input_loc >= 360.0)
			{
				pptr->input_loc -= 360.0;
			}
		}	
		set_input(dptr, pptr);
		
		calculate_activations_excitatory(dptr, pptr, cptr, bptr, timestep);
		calculate_rates_excitatory(dptr, pptr);
		
		
#pragma omp for private (connection)
		for(connection=0; connection<pptr->num_E1_cells; connection++) 
		{
			fill_RC_buffer(dptr, pptr, bptr, connection);
		}
		
#pragma omp single
		{
			if(timestep % 100 == 0)
			{
				record_activations_rates(dptr, pptr, (timestep/100));
				record_input(dptr,pptr,(timestep/100));
				load_full_RC(cptr, pptr);
				log_weight_vector(dptr, pptr, cptr, (timestep/100));
								
			}
		}			
	}
	
	
	
	for(timestep=pptr->start_length; timestep<pptr->training_timesteps; timestep++)
	{
#pragma omp single nowait
	{
		printf("\rTimestep %d", timestep+1);
		fflush(stdout);
	}
		

		record_previous_states(dptr, pptr, cptr);
		
		
#pragma omp single
	{
		increment = pptr->velocity*pptr->timestep_size;
		pptr->input_loc += increment;
		if(pptr->input_loc >= 360.0)
		{
			pptr->input_loc -= 360.0;
		}
	}	
		set_input(dptr, pptr);
		
		calculate_activations_excitatory(dptr, pptr, cptr, bptr, timestep);
		calculate_rates_excitatory(dptr, pptr);
		

#pragma omp for private (connection)
		for(connection=0; connection<pptr->num_E1_cells; connection++) 
		{
			fill_RC_buffer(dptr, pptr, bptr, connection);
		}

		update_RC_weights(dptr, cptr, pptr, bptr);

		
#pragma omp single
		{
			if(timestep % 100 == 0)
			{
				record_activations_rates(dptr, pptr, (timestep/100));
				record_input(dptr,pptr,(timestep/100));
				load_full_RC(cptr, pptr);
				log_weight_vector(dptr, pptr, cptr, (timestep/100));
		
			}
		}			
	}
	
	
#pragma omp single nowait
	{
		printf("\n\nTesting Phase...\n");
		fflush(stdout);
	}
	
		zero_input(dptr, pptr);
	
	if(pptr->symtest)
	{

#pragma omp single
		{
			load_full_RC(cptr, pptr);
			save_RCweights_intermediate(cptr,pptr);
			calculate_weight_vector(dptr, pptr, cptr, 1);

		}
		
		overlay_symmetrical(pptr, dptr, cptr);
	}
	
#pragma omp single
	{	
		pptr->testphase = 1; 
	}
	
	for(timestep=(pptr->training_timesteps); timestep<pptr->timesteps; timestep++)
	{
#pragma omp single nowait
	{
		printf("\rTimestep %d", timestep+1);
		fflush(stdout);
	}
		
		record_previous_states(dptr, pptr, cptr);
		
		calculate_activations_excitatory(dptr, pptr, cptr, bptr, timestep);

		calculate_rates_excitatory(dptr, pptr);
				
#pragma omp for private (connection)
		for(connection=0; connection<pptr->num_E1_cells; connection++)		
		{
			fill_RC_buffer(dptr, pptr, bptr, connection);
		}
		
	
#pragma omp single
		{
			if(timestep % 100 == 0)
			{
				record_activations_rates(dptr, pptr, (timestep/100));
				record_input(dptr,pptr,(timestep/100));
				load_full_RC(cptr, pptr);
				log_weight_vector(dptr, pptr, cptr, (timestep/100));
				
		
			}
		
		
		}		
		
	}
}
	
	end = omp_get_wtime();
	
	printf("Runtime = %f\n", end-start);
	
	
	printf("\n\nSaving results...\n");
	fflush(stdout);

	
	
	save_input(dptr,pptr);
		

	save_rates_activations(dptr, pptr);
	
	load_full_RC(cptr, pptr);
	
	save_RCweights_final(cptr, pptr);
	calculate_weight_vector(dptr, pptr, cptr, 2);

	save_weight_vectors(dptr, pptr);
	
	save_wVecs(dptr, pptr);
	
	free_memory(dptr, cptr, bptr, pptr);
	
	free(dptr);
	free(cptr);
	free(pptr);
	
	for(connection=0;connection<pptr->num_E1_cells; connection++)
	{
		free(bptr[connection]);
	}
	
	free(bptr);
	
	
	printf("\nProgram successfully finished\n\n");
	fflush(stdout);
	
	return 0;
	
}



