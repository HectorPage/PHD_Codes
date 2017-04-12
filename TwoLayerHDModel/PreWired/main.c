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
	int counter;
	int connection;
	double start, end;
	params *pptr;
	data *dptr;
	connect *cptr;
	cArray **combptr;
	cArray **hdptr;

	if(argc!=3)
	{
		fprintf(stderr, "Usage: Oscillate -fname <parameters_file>\n");
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
	
	printf("Parameters read....\n");
	fflush(stdout);
	
	combptr = (cArray **)malloc(sizeof(cArray *) * pptr->num_COMB_cells);
	
	hdptr = (cArray **)malloc(sizeof(cArray *) * pptr->num_HD_cells);
		
	for(counter=0; counter<pptr->num_COMB_cells; counter++) 
	{
		combptr[counter] = malloc(sizeof(cArray));
	}
	
	for(counter=0; counter<pptr->num_HD_cells; counter++) 
	{
		hdptr[counter] = malloc(sizeof(cArray));
	}
		
	build_network(dptr, pptr, cptr, combptr, hdptr);
	
	printf("Network memory allocated....\n");
	fflush(stdout);
	
	set_connectivity(cptr, pptr);
	
	printf("Connectivity set....\n");
	fflush(stdout);
	
    set_favoured_view(dptr, pptr);
	
	printf("Favoured view set....\n");
	fflush(stdout);
	
	set_HD_COMB_weights(pptr,dptr,cptr);
	
	printf("HD to COMB weights set....\n");
	fflush(stdout);

	set_COMB_HD_weights(pptr,dptr,cptr);
	
	printf("COMB to HD weights set....\n");
	fflush(stdout);

	/* Simulating original cue conflict experiment */
	
	zero_states(dptr, pptr, combptr, hdptr);
	
	printf("States zeroed...\n");
	fflush(stdout);

	zero_COMB(dptr, pptr);
	zero_initialiser(dptr, pptr);
	


	printf("Network now being simulated...\n");
	fflush(stdout);
	
	/*This is an initial presentation of the light, to set up HD memory packet*/
	
	printf("INITIALISING HD PACKET\n");
	fflush(stdout);
	
	//PARALLEL REGION BEGINS HERE
	
	omp_set_num_threads(pptr->numThreads);
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
		
#pragma omp single
		{
			dptr->rates_NOROT = 1.0;
			dptr->rates_ROT = 0.0;
		}
		
		set_initialiser(dptr, pptr);
				
		calculate_activations_excitatory(dptr, pptr, cptr, combptr, hdptr, timestep);
		calculate_rates_excitatory(dptr, pptr);	
		
#pragma omp for private(connection)
		for(connection=0; connection<pptr->num_COMB_cells; connection++)
		{
			fill_comb_buffer(dptr, pptr, combptr, connection);
		}
		
#pragma omp for private(connection)
		for(connection=0; connection<pptr->num_HD_cells; connection++)
		{
			fill_hd_buffer(dptr, pptr, hdptr, connection);
		}
		
#pragma omp single
		{
			if(timestep % 100 == 0)
			{
				record_activations_rates(dptr, pptr, (timestep/100));	//If uncommented it is timestep/100
				record_initialiser(dptr,pptr,(timestep/100));
			}
		}			
	}
	
	zero_initialiser(dptr, pptr);
		
#pragma	omp single nowait
		{
			printf("\nIntermediate ff_weights saved....\n");
			fflush(stdout);
		}
		
	calculate_pvector_initialiser(dptr, pptr); 
	
#pragma	omp single nowait
		{
			printf("Conflict pvector calculated......\n");
			fflush(stdout);
		}
		
#pragma omp single
		{
			save_pvector(dptr);
		}
	
	/*HOLDING PACKET*/
	
#pragma	omp single nowait
		{
			printf("Holding activity packet......\n");
			fflush(stdout);
		}
		
	for(timestep = pptr->input_timesteps; timestep<(pptr->input_timesteps+pptr->pausetimesteps); timestep++)
	{
#pragma	omp single nowait
		{
			printf("\rTimestep %d", timestep+1);
			fflush(stdout);
		}
		
		record_previous_states(dptr, pptr, cptr);
		
#pragma omp single
		{
			dptr->rates_NOROT = 1.0;
			dptr->rates_ROT = 0.0;
		}
		
		calculate_activations_excitatory(dptr, pptr, cptr, combptr, hdptr, timestep);
		calculate_rates_excitatory(dptr, pptr);	
		
#pragma omp for private(connection)
		for(connection=0; connection<pptr->num_COMB_cells; connection++)
		{
			fill_comb_buffer(dptr, pptr, combptr, connection);
		}
		
#pragma omp for private(connection)
		for(connection=0; connection<pptr->num_HD_cells; connection++)
		{
			fill_hd_buffer(dptr, pptr, hdptr, connection);
		}
		
		
#pragma omp single
		{
			if(timestep % 100 == 0)
			{
				record_activations_rates(dptr, pptr, (timestep/100));
				record_initialiser(dptr,pptr,(timestep/100));
			}
		}
		
	}
	

	
	/*PATH INTEGRATION*/
#pragma	omp single nowait
		{
			printf("\nROTATING RAT\n");
			fflush(stdout);
		}
		
	for(timestep = (pptr->input_timesteps + pptr->pausetimesteps) ; timestep<(pptr->rotation_timesteps + pptr->input_timesteps + pptr->pausetimesteps); timestep++)
	{
#pragma	omp single nowait
		{
			printf("\rTimestep %d", timestep+1);
			fflush(stdout);
		}
		
		record_previous_states(dptr, pptr, cptr);
		
#pragma omp single
		{
			dptr->rates_ROT = 1.0;
			dptr->rates_NOROT = 0.0;
		}
		
		calculate_activations_excitatory(dptr, pptr, cptr, combptr, hdptr, timestep);
		calculate_rates_excitatory(dptr, pptr);
		
#pragma omp for private(connection)
		for(connection=0; connection<pptr->num_COMB_cells; connection++)
		{
			fill_comb_buffer(dptr, pptr, combptr, connection);
		}
#pragma omp for private(connection)
		for(connection=0; connection<pptr->num_HD_cells; connection++)
		{
			fill_hd_buffer(dptr, pptr, hdptr, connection);
		}
		
		
#pragma omp single
		{
			if(timestep % 100 == 0)
			{
				record_activations_rates(dptr, pptr, (timestep/100));
				record_initialiser(dptr,pptr,(timestep/100));
			
			}
		}	
	}
	
	
	/*HOLDING PACKET*/
#pragma	omp single nowait
		{
			printf("\nHolding activity packet\n");
			fflush(stdout);
		}
		
	for(timestep=(pptr->rotation_timesteps + pptr->input_timesteps + pptr->pausetimesteps); timestep<pptr->timesteps; timestep++)
	{
#pragma	omp single nowait
		{
			printf("\rTimestep %d", timestep+1);
			fflush(stdout);
		}
		
#pragma omp single
		{
			dptr->rates_ROT = 0.0;
			dptr->rates_NOROT = 1.0;
		}
		
		record_previous_states(dptr, pptr, cptr);
		
		
		calculate_activations_excitatory(dptr, pptr, cptr, combptr, hdptr, timestep);
		calculate_rates_excitatory(dptr, pptr);
		
#pragma omp for private(connection)
		for(connection=0; connection<pptr->num_COMB_cells; connection++)
		{
			fill_comb_buffer(dptr, pptr, combptr, connection);
		}
#pragma omp for private(connection)
		for(connection=0; connection<pptr->num_HD_cells; connection++)
		{
			fill_hd_buffer(dptr, pptr, hdptr, connection);
		}
		
		
		
#pragma omp single
		{
			if(timestep % 100 == 0)
			{
				record_activations_rates(dptr, pptr, (timestep/100));
				record_initialiser(dptr,pptr,(timestep/100));
			
			}
		}
		
	}

	
}
	
	end = omp_get_wtime();
	printf("\nRuntime = %f\n", end-start);
	
	/*SAVING RESULTS*/
		
	
	printf("\nSAVING RESULTS\n");
	fflush(stdout);
	
	save_initialiser(dptr,pptr);
	
	printf("Saved HD packet locations\n");
	fflush(stdout);
	
	save_rates_activations(dptr, pptr);
	
	printf("Saved rates and activations\n");
	
	free_memory(pptr, dptr, cptr, combptr, hdptr);
	free(dptr);
	free(cptr);
	free(pptr);
	
	for(connection=0;connection<pptr->num_COMB_cells; connection++)
	{
		free(combptr[connection]);
	}
	
	free(combptr);
	
	
	for(connection=0;connection<pptr->num_HD_cells; connection++)
	{
		free(hdptr[connection]);
	}	
	
	free(hdptr);
	
	printf("Program successfully finished\n");
	fflush(stdout);
	
	return 0;
	
}



