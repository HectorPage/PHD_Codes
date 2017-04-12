#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
	int idx;
	
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
	
	read_parameters(fname, pptr, cptr);
	
	dptr = (data *)malloc(sizeof(data));
	
	build_network(dptr, pptr, cptr);
	
	set_connectivity(cptr, pptr);
	
    set_favoured_view(dptr, pptr);
	
	set_RC_weights(pptr, dptr, cptr);
	
	for(idx=0; idx<10; idx++)
	{
		printf("%f\n", cptr->weights_visring[0][idx]);
	}
	
	set_visring_weights(pptr,dptr,cptr);
	
	for(idx=0; idx<10; idx++)
	{
		printf("%f\n", cptr->weights_visring[0][idx]);
	}
	
	
	save_ffweights_init(cptr, pptr);
	
	/* Simulate Network */
	
	zero_states(dptr, pptr);
	zero_vis1(dptr, pptr);
	zero_visring(dptr,pptr);
	zero_averages(dptr, pptr);

		
	printf("\nNetwork now being simulated...\n\n");
	fflush(stdout);
	
	for(timestep=0; timestep<pptr->input_timesteps; timestep++)
	{
		printf("\rTimestep %d", timestep+1);
		fflush(stdout);
		
		record_previous_states(dptr, pptr, cptr);
				
		set_vis1(dptr, pptr);
				
		calculate_activations_excitatory(dptr, pptr, cptr, timestep);
		calculate_rates_excitatory(dptr, pptr);
		calculate_E1_rates_average(dptr, pptr);
		
		if(timestep % 100 == 0)
		{
		record_activations_rates(dptr, pptr, (timestep/100));
		record_vis1(dptr,pptr,(timestep/100));
		record_visring(dptr,pptr,(timestep/100));
		}
				
	}
	
	zero_vis1(dptr, pptr);
	
	for(timestep=pptr->input_timesteps; timestep<pptr->timesteps; timestep++)
	{
		printf("\rTimestep %d", timestep+1);
		fflush(stdout);
		
		record_previous_states(dptr, pptr, cptr);
		
		set_visring(dptr, pptr);
		
		calculate_activations_excitatory(dptr, pptr, cptr, timestep);
		calculate_rates_excitatory(dptr, pptr);
		calculate_E1_rates_average(dptr, pptr);
		
		if(pptr->learning)
		{
		update_visring_weights(dptr, cptr, pptr);
		}
		
		if(timestep % 100 == 0)
		{
		record_activations_rates(dptr, pptr, (timestep/100));
		record_vis1(dptr,pptr,(timestep/100));
		record_visring(dptr,pptr,(timestep/100));
		}
		
		
	}		
		
	
	printf("\n\nSaving results...\n");
	fflush(stdout);
	
	calculate_pvector_vis1(dptr, pptr);
	
	save_pvector(dptr);
	save_visinput1(dptr,pptr);	
	save_visring(dptr,pptr);
	
	save_rates_activations(dptr, pptr);
	
	save_ffweights_final(cptr, pptr);
	
	free_memory(dptr, cptr);
	free(dptr);
	free(cptr);
	free(pptr);
	
	printf("\nProgram successfully finished\n\n");
	fflush(stdout);
	
	return 0;
	
}



