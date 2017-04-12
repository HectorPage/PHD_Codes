#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "calculate.h"
#include "data.h"
#include "file_io.h"
#include "init.h"
#include <omp.h>


int main(int argc, char **argv)
{
	char *fname;
	int timestep;
	float increment;
	double start, end;
	params *pptr;
	data *dptr;
	connect *cptr;
	if(argc!=3)
	{
		fprintf(stderr, "Usage: nopi -fname <parameters_file>\n");
		exit(1);
	}
	else
	{
		pptr = (params *)malloc(sizeof(params));
		cptr = (connect *)malloc(sizeof(connect));
		fname = argv[2];
	}
	
	/*Setting up network*/
	
	dptr = (data *)malloc(sizeof(data));
	
	read_parameters(fname, pptr, cptr);
	
	printf("\nParameters read....\n\n");
	fflush(stdout);
	
	build_network(dptr, pptr, cptr);
	
	printf("\nNetwork memory allocated....\n\n");
	fflush(stdout);
	
	set_connectivity(cptr, pptr);
	
    set_favoured_view(dptr, pptr);
	
	set_RC_weights(pptr, dptr, cptr);
	
	set_visring_weights(pptr,dptr,cptr);
	
	save_ffweights_init(cptr, pptr);
	
	calculate_weight_vector(dptr, pptr, cptr, 0);
		
	printf("\nWeights initialised....\n\n");
	fflush(stdout);
	
	/* Simulating original cue conflict experiment */
	
	zero_states(dptr, pptr);
	zero_pi(dptr, pptr);
	zero_visring(dptr,pptr);

	printf("\nNetwork now being simulated...\n\n");
	fflush(stdout);
	
	/*INITIAL CUE CONFLICT*/
	
	/*This is an initial presentation of the light, to set up HD memory packet*/
	
	printf("\nINITIALISING HD PACKET\n\n");
	fflush(stdout);
	
	
	omp_set_num_threads(pptr->threads);
	
	start = omp_get_wtime();
	
#pragma omp parallel private (timestep)
{
	
	for(timestep=0; timestep<pptr->input_timesteps; timestep++)
	{
#pragma omp single nowait
		{
			printf("\rTimestep %d", timestep+1);
			fflush(stdout);
		}
		
		record_previous_states(dptr, pptr, cptr);
	
		initialise_pi(dptr, pptr);
				
		calculate_activations_excitatory(dptr, pptr, cptr, timestep);
		calculate_rates_excitatory(dptr, pptr);
		
#pragma omp single
		{
			if(timestep % 100 == 0)
			{
				record_activations_rates(dptr, pptr, (timestep/100));
				record_pi(dptr,pptr,(timestep/100));
				record_visring(dptr,pptr,(timestep/100));
			}
		}		
	}
	
	zero_pi(dptr, pptr);
	
	/*This is now generating a cue conflict situation*/
	
#pragma omp single nowait
	{
		printf("\nGENERATING CUE CONFLICT\n\n");
		fflush(stdout);
	}
	
	for(timestep=pptr->input_timesteps; timestep<pptr->conflict_timesteps; timestep++)
	{
#pragma omp single nowait
		{
			printf("\rTimestep %d", timestep+1);
			fflush(stdout);
		}
		
		record_previous_states(dptr, pptr, cptr);
		
		set_visring(dptr, pptr);
		
		calculate_activations_excitatory(dptr, pptr, cptr, timestep);
		calculate_rates_excitatory(dptr, pptr);
				
		update_visring_weights(dptr, cptr, pptr);
		
#pragma omp single
		{
			if(timestep % 100 == 0)
			{
				record_activations_rates(dptr, pptr, (timestep/100));
				record_pi(dptr,pptr,(timestep/100));
				record_visring(dptr,pptr,(timestep/100));
			}
		}
		
	}	

#pragma omp single
{
	save_ffweights_intermediate(cptr, pptr); /*Saves feed-forward weights after the initial cue conflict*/
	
	printf("\nIntermediate ff_weights saved....\n\n");
	fflush(stdout);
		
	calculate_pvector(dptr, pptr); 
	
	printf("\nConflict pvector calculated......\n\n");
	fflush(stdout);
	
	save_pvector(dptr);
	
	printf("\nConflict pvector saved....\n\n");
	fflush(stdout);
	
	pptr->pi_loc = dptr->pvector;
	
	calculate_weight_vector(dptr, pptr, cptr, 1);
		
	printf("\nIntermediate weight vector calculated....\n\n");
	fflush(stdout);

	
	/*PATH INTEGRATION*/
	
	printf("\n\nROTATING RAT\n\n");
	fflush(stdout);
}	
	for(timestep = pptr->conflict_timesteps; timestep<pptr->pi_timesteps; timestep++)
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
		pptr->visring_loc += increment;
		if(pptr->visring_loc >= 360.0)
		{
			pptr->visring_loc -= 360.0;
		}
  	}
		set_visring(dptr, pptr); 
		
		
		/*FORCING HD TO MOVE AT CORRECT SPEED USING INITIALISER*/
	
#pragma omp single
		{
			
			pptr->pi_loc += increment;
			if(pptr->pi_loc >= 360.0)
			{
				pptr->pi_loc -= 360.0;
			}
		}
		
			set_pi(dptr, pptr);
		
		
		calculate_activations_excitatory(dptr, pptr, cptr, timestep);
		calculate_rates_excitatory(dptr, pptr);
				
		
		update_visring_weights(dptr, cptr, pptr);
		
#pragma omp single
		{
			if(timestep % 100 == 0)
			{
				record_activations_rates(dptr, pptr, (timestep/100));
				record_pi(dptr,pptr,(timestep/100));
				record_visring(dptr,pptr,(timestep/100));
			}
		}
				
	}
	
	
	/*Attractor section*/
#pragma omp single nowait
	{
		printf("\n\nHOLDING PACKET\n");
		fflush(stdout);
	}
	
	zero_pi(dptr, pptr);
	zero_visring(dptr, pptr);
	
	for(timestep=pptr->pi_timesteps;timestep<(pptr->pi_timesteps+10000); timestep++)
	{
#pragma omp single nowait
		{
			printf("\rTimestep %d", timestep+1);
			fflush(stdout);
		}
		
		record_previous_states(dptr, pptr, cptr);
		
		calculate_activations_excitatory(dptr, pptr, cptr, timestep);
		calculate_rates_excitatory(dptr, pptr);
		
#pragma omp single
		{
			if(timestep % 100 == 0)
			{
				record_activations_rates(dptr, pptr, (timestep/100));
				record_pi(dptr,pptr,(timestep/100));
				record_visring(dptr,pptr,(timestep/100));
			}
		}
		
		
	}
	
	/*Extra test section - following the light cue*/
	
#pragma omp single
	{
		printf("\n\nTRACKING LIGHT\n");
		fflush(stdout);
	}
	
	calculate_pvector(dptr,pptr);
	pptr->visring_loc = dptr->pvector;
	
	for(timestep=(pptr->pi_timesteps+10000);timestep<pptr->timesteps; timestep++)
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
		
			pptr->visring_loc += increment;
			if(pptr->visring_loc >= 360.0)
			{
				pptr->visring_loc -= 360.0;
			}
		}
		
		set_visring(dptr, pptr);
		
		calculate_activations_excitatory(dptr, pptr, cptr, timestep);
		calculate_rates_excitatory(dptr, pptr);
		
#pragma omp single
		{
			if(timestep % 100 == 0)
			{
				record_activations_rates(dptr, pptr, (timestep/100));
				record_pi(dptr,pptr,(timestep/100));
				record_visring(dptr,pptr,(timestep/100));
			}
		}
		
	}
}
	end = omp_get_wtime();
	
	printf("Runtime = %f\n", end-start);
	
	
	/*SAVING RESULTS*/
		
	
	printf("\n\nSAVING RESULTS\n");
	fflush(stdout);
	
	save_pi_input(dptr,pptr);
	
	printf("\n\nSaved HD packet locations\n\n");
	fflush(stdout);

	save_visring(dptr,pptr);
	
	printf("\n\nSaved visring packet locations\n\n");
	fflush(stdout);

	
	save_rates_activations(dptr, pptr);
	
	printf("\n\nSaved rates and activations\n\n");
	
	save_ffweights_final(cptr, pptr);
	
	printf("\n\nSaved final ff weights\n\n");
	
	calculate_weight_vector(dptr, pptr, cptr, 2);
	
	printf("\n\nCell 350 final weight vector: %f\n\n", dptr->HD_weight_vectors_final[350]);
	
	save_weight_vector(dptr, pptr);
	
	printf("\n\nSaved final ff weight vector\n\n");
	
	free_memory(pptr, dptr, cptr);
	free(dptr);
	free(cptr);
	free(pptr);
	
	printf("\nProgram successfully finished\n\n");
	fflush(stdout);
	
	return 0;
	
}



