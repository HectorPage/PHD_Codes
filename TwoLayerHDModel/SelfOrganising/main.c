#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "calculate.h"
#include "data.h"
#include "file_io.h"
#include "init.h"
#include "array_utils.h"
#include "utils.h"



int main(int argc, char **argv)
{
	char *fname;
	int timestep;
	int counter;
	int epochs;
	int connection;
	int cell;
	float increment;
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
	
	set_ROT_COMB_weights(pptr,dptr,cptr);
	
	printf("ROT to COMB weights set....\n");
	fflush(stdout);
	
	load_full_HDtoCOMB(cptr,pptr);
	save_initial_HDtoCOMB_weights(cptr, pptr);			
	save_initial_COMBtoHD_weights(cptr, pptr);
	save_initial_ROTtoCOMB_weights(cptr, pptr);
	
	/* Simulating original cue conflict experiment */
	
	zero_states(dptr, pptr, combptr, hdptr);
	
	printf("States zeroed...\n");
	fflush(stdout);

	zero_COMB(dptr, pptr);
	zero_initialiser(dptr, pptr);
	
	zero_previous_states(pptr, cptr, dptr);			//New function I added to help debug without needed to record_previous_states

	printf("Network now being simulated...\n");
	fflush(stdout);
	
	//PARALLEL REGION BEGINS HERE
	
	omp_set_num_threads(pptr->numThreads);
	start = omp_get_wtime();
	
	
	//TRAINING PHASE
	
	
	//THIS BIT TRAINS NON-OFFSET CONNECTIVITY

		
		//ThreadNum = omp_get_thread_num();		
//#pragma	omp single
//		{
			printf("\nTRAINING PHASE....\n");
			fflush(stdout);
//		}

		//here loop over epochs

		for(epochs=1;epochs<(pptr->training_loops+1);epochs++) 
		{
			printf("\nEpoch %d \n", epochs);
			fflush(stdout);

			
#pragma omp parallel private (timestep) //keep epochs private
			{
#pragma omp single
			{
			cell=1;
			}
			//here loop over timesteps (dependent on epoch), with timestep private

		for(timestep=(((epochs-1)*(pptr->NonOffTrainTimesteps*pptr->num_HD_cells)) + ((epochs-1)*pptr->rotation_timesteps));timestep<((epochs*(pptr->NonOffTrainTimesteps*pptr->num_HD_cells)) + ((epochs-1)*pptr->rotation_timesteps));timestep++)
		{
					
/*#pragma omp single nowait
//if(ThreadNum==0)
			{
				printf("\rTimestep %d", timestep+1);
				fflush(stdout);
			}
	*/		
#pragma omp single			
//if(ThreadNum==0)
{
	if(pptr->introtest)
		{

			
				if(timestep == pptr->conduction_buffer_size +1)
				{
					FILE *fptr;
					fptr = file_open("_COMBslice.bdat", "wb");
					
					fwrite(&(dptr->rates_COMB[0]), sizeof(float), pptr->num_COMB_cells, fptr);
					
					fclose(fptr);
					
					
					fptr = file_open("_COMBactslice.bdat", "wb");
					
					fwrite(&(dptr->activations_COMB[0]), sizeof(float), pptr->num_COMB_cells, fptr);
					
					fclose(fptr);
					
					
					
					
					fptr = file_open("_HDslice.bdat", "wb");
					
					fwrite(&(dptr->rates_HD[0]), sizeof(float), pptr->num_HD_cells, fptr);
					
					fclose(fptr);
					
					
					exit(0);
				}
				
			}
}

			
			
		
			
			pptr->initialiser_loc = dptr->favoured_view[cell];
			set_initialiser(dptr, pptr);
			
#pragma omp for private(connection)			
			for(connection=0;connection<pptr->num_ROT_cells/2;connection++)
				{
					dptr->rates_ROT[connection] = 1.0;
				}
			
#pragma omp for private(connection)				
			for(connection=pptr->num_ROT_cells/2;connection<pptr->num_ROT_cells;connection++)
				{
					dptr->rates_ROT[connection] = 0.0;
				}
					
			
			
			record_previous_states(dptr, pptr, cptr, hdptr);
				
			
			calculate_activations_excitatory(dptr, pptr, cptr, combptr, hdptr);	
			
			
			calculate_rates_excitatory(dptr, pptr);	

#pragma omp single
//if(ThreadNum==0)
			{
				if(timestep % 100 == 0)
				{
					record_activations_rates(dptr, pptr, (timestep/100));
					record_initialiser(dptr,pptr,(timestep/100));
				}
			}	

			
			
			update_HDtoCOMB_weights(dptr, cptr, pptr, hdptr);
			update_COMBtoHD_weights(dptr, cptr, pptr, combptr);		
			
			
			update_ROTtoCOMB_weights(dptr, cptr, pptr);
					
						
#pragma omp single
//if(ThreadNum==0)
{
				
				if(timestep%pptr->NonOffTrainTimesteps == 0 && timestep>0)
				{
					cell+= 1;
					if (cell>= pptr->num_HD_cells)
						cell -= (pptr->num_HD_cells - 1);
				}
				
			}
			
			
			
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
			

	}
		
		
		
	//THIS BIT TRAINS OFFSET CONNECTIVITY
	

			
	for(timestep=((epochs*(pptr->NonOffTrainTimesteps*pptr->num_HD_cells)) + ((epochs-1)*pptr->rotation_timesteps));timestep<((epochs*(pptr->NonOffTrainTimesteps*pptr->num_HD_cells)) + (epochs*pptr->rotation_timesteps));timestep++)
	{
#pragma omp for private(connection)				
		for(connection=0;connection<pptr->num_ROT_cells/2;connection++)
		{
			dptr->rates_ROT[connection] = 0.0;
		}
		
#pragma omp for private(connection)		
		for(connection=pptr->num_ROT_cells/2;connection<pptr->num_ROT_cells;connection++)
		{
			dptr->rates_ROT[connection] = 1.0;
		}
		
		
/*#pragma	omp single nowait
//if(ThreadNum==0)
{
			printf("\rTimestep %d", timestep+1);
			fflush(stdout);
		}
*/		
		
	
		
#pragma omp single
//if(ThreadNum==0)
{
			increment = pptr->velocity*pptr->timestep_size;
			pptr->initialiser_loc += increment;
			
			if(pptr->initialiser_loc >= 360.0)
			{
				pptr->initialiser_loc -= 360.0;
			}
		}
		
		
		set_initialiser(dptr, pptr);	
		
		record_previous_states(dptr, pptr, cptr, hdptr);
		
		calculate_activations_excitatory(dptr, pptr, cptr, combptr, hdptr);
		
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
//if(ThreadNum==0)
{
			if(timestep % 100 == 0)
			{
				record_activations_rates(dptr, pptr, (timestep/100));
				record_initialiser(dptr,pptr,(timestep/100));
			}
		}				
		
		
		update_HDtoCOMB_weights(dptr, cptr, pptr, hdptr);	
		
		update_COMBtoHD_weights(dptr, cptr, pptr, combptr);
		
		update_ROTtoCOMB_weights(dptr, cptr, pptr);
		

			}
		}	
		}
	
#pragma omp parallel private(timestep)
	{
	
//NOW THE TESTING PHASE
		
#pragma omp single
		{
			load_full_HDtoCOMB(cptr,pptr);
			save_HDtoCOMB_weights(cptr, pptr);			
			save_COMBtoHD_weights(cptr, pptr);
			save_ROTtoCOMB_weights(cptr, pptr);
	
		}
		
#pragma	omp single
		{
			printf("\nTESTING PHASE....\n");
			fflush(stdout);
		}
	
	
	zero_initialiser(dptr, pptr);
	pptr->extern_inhib = 0.0;	//This is the feed-forward inhibition.

	
	/*HOLDING PACKET*/
	
#pragma	omp single
		{
			printf("Holding activity packet......\n");
			fflush(stdout);
		}
		
		
	for(timestep = ((pptr->NonOffTrainTimesteps*pptr->num_HD_cells*pptr->training_loops) + (pptr->training_loops*pptr->rotation_timesteps)); timestep<(((pptr->NonOffTrainTimesteps*pptr->num_HD_cells*pptr->training_loops) + (pptr->rotation_timesteps*pptr->training_loops))+pptr->pausetimesteps); timestep++)
	{
/*#pragma	omp single nowait
		{
			printf("\rTimestep %d", timestep+1);
			fflush(stdout);
		}
*/		
#pragma omp for private(connection)		
		for(connection=0;connection<pptr->num_ROT_cells/2;connection++)
		{
			dptr->rates_ROT[connection] = 1.0;
		}
		
#pragma omp for private(connection)		
		for(connection=pptr->num_ROT_cells/2;connection<pptr->num_ROT_cells;connection++)
		{
			dptr->rates_ROT[connection] = 0.0;
		}
		
		record_previous_states(dptr, pptr, cptr, hdptr);
		
		
		calculate_activations_excitatory(dptr, pptr, cptr, combptr, hdptr);
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
#pragma	omp single
		{
			printf("\nROTATING RAT\n");
			fflush(stdout);
		}
		
	for(timestep = ((pptr->NonOffTrainTimesteps*pptr->num_HD_cells*pptr->training_loops) + (pptr->rotation_timesteps*pptr->training_loops)) + pptr->pausetimesteps ; timestep<(((pptr->NonOffTrainTimesteps*pptr->num_HD_cells*pptr->training_loops) + (pptr->rotation_timesteps*pptr->training_loops))+pptr->pausetimesteps) + pptr->test_timesteps; timestep++)
	{
/*#pragma	omp single nowait
		{
			printf("\rTimestep %d", timestep+1);
			fflush(stdout);
		}
*/		
#pragma omp for private(connection)		
		for(connection=0;connection<pptr->num_ROT_cells/2;connection++)
		{
			dptr->rates_ROT[connection] = 0.0;
		}
		
#pragma omp for private(connection)		
		for(connection=pptr->num_ROT_cells/2;connection<pptr->num_ROT_cells;connection++)
		{
			dptr->rates_ROT[connection] = 1.0;
		}		
		record_previous_states(dptr, pptr, cptr, hdptr);
		
		
		calculate_activations_excitatory(dptr, pptr, cptr, combptr, hdptr);
		
		
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
#pragma	omp single
		{
			printf("\nHolding activity packet\n");
			fflush(stdout);
		}
		
	for(timestep=(((pptr->NonOffTrainTimesteps*pptr->num_HD_cells*pptr->training_loops) + (pptr->rotation_timesteps*pptr->training_loops))+pptr->pausetimesteps) + pptr->test_timesteps; timestep<pptr->timesteps; timestep++)
	{
/*#pragma	omp single nowait
		{
			printf("\rTimestep %d", timestep+1);
			fflush(stdout);
		}
*/
 #pragma omp for private(connection)
		for(connection=0;connection<pptr->num_ROT_cells/2;connection++)
		{
			dptr->rates_ROT[connection] = 1.0;
		}
		
#pragma omp for private(connection)		
		for(connection=pptr->num_ROT_cells/2;connection<pptr->num_ROT_cells;connection++)
		{
			dptr->rates_ROT[connection] = 0.0;
		}	
		
		record_previous_states(dptr, pptr, cptr, hdptr);
		
		
		
		calculate_activations_excitatory(dptr, pptr, cptr, combptr, hdptr);
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
	
	printf("Memory freed\n");

	
	printf("Program successfully finished\n");
	fflush(stdout);
	
	return 0;
	
}



