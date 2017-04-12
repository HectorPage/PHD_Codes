#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "calculate.h"
#include "data.h"

void zero_states(data *dptr, params *pptr, cArray **bptr, cArray **bxptr)
{
	int cell;
	int connection;
	
	for(cell=0; cell<pptr->num_HD_cells; cell++)
	{
		dptr->activations_HD[cell] = 0.0;
		dptr->prev_activations_HD[cell] = 0.0;
		dptr->rates_HD[cell] = 0.0;
		dptr->prev_rates_HD[cell] = 0.0;
	}
	
	for(cell=0; cell<pptr->conduction_buffer_size; cell++)
	{
		for(connection=0; connection<pptr->num_COMB_cells; connection++)
		{
			
		bptr[connection]->array[cell] = 0.0;
				
		}
	}
	
	printf("COMB buffer cleared\n");
	fflush(stdout);

	
	for(cell=0; cell<pptr->conduction_buffer_size; cell++)
	{
		for(connection=0; connection<pptr->num_HD_cells; connection++)
		{
			
			bxptr[connection]->array[cell] = 0.0;
			
		}
	}
	
	//NOT SURE ABOUT THE ABOVE ZEROING OF BUFFERS - MAY NEED TO SWAP NUM_HD AND NUM_COMB?
	
	
	dptr->rates_ROT = 0.0;
	dptr->rates_NOROT = 0.0;

	
		return;
}

void zero_COMB(data *dptr, params *pptr)
{
	int cell;
	
	
	for(cell=0; cell<pptr->num_COMB_cells; cell++)
	{
		dptr->activations_COMB[cell] = 0.0;
		dptr->prev_activations_COMB[cell] = 0.0;
		dptr->rates_COMB[cell] = 0.0;
		dptr->prev_rates_COMB[cell] = 0.0;
	}
	
	
	return;
}


void set_initialiser(data *dptr, params *pptr)
{
    int cell;
    float distance, distance1, distance2;
    
#pragma omp for private(cell, distance, distance1, distance2)
		for(cell=0; cell<pptr->num_HD_cells; cell++)
		{
			distance1 = fabs(pptr->initialiser_loc - dptr->favoured_view[cell]);
			distance2 = fabs(360.0 - distance1);
        
			if(distance1<=distance2)
				distance = distance1;
			else
				distance = distance2;
        
			dptr->initialiser[cell] = pptr->initialiser_str * exp(-0.5 * (distance/pptr->initialiser_sigma) * (distance/pptr->initialiser_sigma));
		}
	
}


void zero_initialiser(data *dptr, params *pptr)
{
	int cell;
#pragma omp for private(cell)
	for(cell=0; cell<pptr->num_HD_cells; cell++)
	{
		dptr->initialiser[cell] = 0.0;
	}
	
	return;
}


void record_previous_states(data *dptr, params *pptr, connect *cptr)
{
#pragma omp sections
	{
#pragma omp section
		{
			memcpy(&(dptr->prev_activations_HD[0]), &(dptr->activations_HD[0]), pptr->num_HD_cells*sizeof(float));
		}
#pragma omp section
		{	
			memcpy(&(dptr->prev_rates_HD[0]), &(dptr->rates_HD[0]), pptr->num_HD_cells*sizeof(float));
		}
#pragma omp section
		{
			memcpy(&(dptr->prev_activations_COMB[0]), &(dptr->activations_COMB[0]), pptr->num_COMB_cells*sizeof(float));
		}
#pragma omp section
		{
			memcpy(&(dptr->prev_rates_COMB[0]), &(dptr->rates_COMB[0]), pptr->num_COMB_cells*sizeof(float));
		}
	}	
	return;
}



void calculate_activations_excitatory(data *dptr, params *pptr, connect *cptr, cArray **combptr, cArray **hdptr, int timestep)
{
	int post_cell, synapse, connection;
	float HD_coefficient, COMB_coefficient;
	float HD_inhibitory, COMB_inhibitory, ROT_input, NOROT_input, COMB_input, HD_input_COMB;
	float HD_COMB_scale, COMB_HD_scale, ROT_COMB_scale;
	
	HD_coefficient = pptr->timestep_size/pptr->tau_HD;
	COMB_coefficient = pptr->timestep_size/pptr->tau_COMB;
		
	HD_COMB_scale = pptr->phi_HDCOMB/(float)pptr->num_HD_cells;
	COMB_HD_scale = pptr->phi_COMBHD/(float)pptr->num_COMB_cells;
	ROT_COMB_scale = pptr->phi_ROT/(float)pptr->num_ROT_cells; //THIS IS ALSO USED FOR NOROT, SINCE THEY ARE PHYSIOLOGICALLY FROM THE SAME SOURCE
	
	/*Calculating Activation for HD cells*/
	
	HD_inhibitory = 0.0;
	
	for(synapse = 0; synapse<pptr->num_HD_cells; synapse++)
	{
		HD_inhibitory += pptr->global_inhibition * dptr->prev_rates_HD[synapse];
	}
#pragma omp for private(post_cell, COMB_input, connection)
	for(post_cell=0; post_cell<pptr->num_HD_cells; post_cell++)
	{
	
		COMB_input = 0.0;
		
		for(connection=0; connection<pptr->num_COMB_cells; connection++)
		{
			COMB_input += cptr->COMB_HD_weights[post_cell][connection] * read_comb_buffer(combptr, connection);
		}
		
				
		dptr->activations_HD[post_cell] = ((1.0 - HD_coefficient) * dptr->prev_activations_HD[post_cell])
		+ (HD_coefficient * dptr->initialiser[post_cell])
		+ (HD_coefficient * (COMB_input * COMB_HD_scale))
		- (HD_coefficient * HD_inhibitory)
		- (HD_coefficient * pptr->extern_inhib);
		

	}
	
	/*Calculating Activation for COMB cells*/
	
		
	COMB_inhibitory = 0.0;
	
	for(synapse=0; synapse<pptr->num_COMB_cells; synapse++)
	{
		COMB_inhibitory += pptr->COMB_inhibition * dptr->prev_rates_COMB[synapse];
	}
	
#pragma omp for private(post_cell, HD_input_COMB, connection)
	for(post_cell=0; post_cell<pptr->num_COMB_cells; post_cell++)
	{
		HD_input_COMB= 0.0;
		
		for(connection=0; connection<pptr->num_HD_cells; connection++)
		{
			HD_input_COMB += cptr->HD_COMB_weights[post_cell][connection] * read_hd_buffer(hdptr, connection);
		}
		
	
		if(post_cell<pptr->num_COMB_cells/2)
		{
			ROT_input = dptr->rates_ROT;
			NOROT_input = 0.0;
			
		}
		else{
			ROT_input = 0.0;
			NOROT_input = dptr->rates_NOROT;
		}
	

		dptr->activations_COMB[post_cell] = ((1.0 - COMB_coefficient) * dptr->prev_activations_COMB[post_cell])
		+ ((HD_input_COMB * HD_COMB_scale) * COMB_coefficient)								
		+ ((ROT_input * ROT_COMB_scale) * COMB_coefficient)
		+ ((NOROT_input * ROT_COMB_scale) * COMB_coefficient)
		- (COMB_inhibitory * COMB_coefficient);
		
	}
		
	return;
}

void calculate_rates_excitatory(data *dptr, params *pptr)
{
	int cell;
#pragma omp for private(cell)
	for(cell=0; cell<pptr->num_HD_cells; cell++)
	{
		dptr->rates_HD[cell] = 1.0 / (1.0 + exp(-2.0 * pptr->beta_HD * (dptr->activations_HD[cell] - pptr->alpha_HD)));
	}
#pragma omp for private(cell)
	for(cell=0; cell<pptr->num_COMB_cells; cell++)
	{	
		dptr->rates_COMB[cell] = 1.0 / (1.0 + exp(-2.0 * pptr->beta_COMB * (dptr->activations_COMB[cell] - pptr->alpha_COMB)));
		
	}
		
	return;
}

void record_activations_rates(data *dptr, params *pptr, int index)
{
	int cell;
	
	for(cell=0; cell<pptr->num_HD_cells; cell++)
	{
		dptr->activations_HD_time[cell][index] = dptr->activations_HD[cell];
		dptr->rates_HD_time[cell][index] = dptr->rates_HD[cell];
	}
	
	for(cell=0; cell<pptr->num_COMB_cells; cell++)
	{
		dptr->activations_COMB_time[cell][index] = dptr->activations_COMB[cell];
		dptr->rates_COMB_time[cell][index] = dptr->rates_COMB[cell];
	}
	
	dptr->rates_ROT_time[index] = dptr->rates_ROT;
	dptr->rates_NOROT_time[index] = dptr->rates_NOROT;

	
	return;
}

void calculate_pvector_initialiser(data *dptr, params *pptr)
{
	
	int cell;
	float vector1 = 0.0;
	float vector2 = 0.0;
	
#pragma omp for private(cell, vector1, vector2)
	for(cell=0; cell<pptr->num_HD_cells; cell++)
	{
		vector1+= dptr->rates_HD[cell] * (sinf(dptr->favoured_view[cell]*(PI/180.0)));
		vector2+= dptr->rates_HD[cell] * (cosf(dptr->favoured_view[cell]*(PI/180.0)));
	}
	
	if (vector1 > 0.0 && vector2 > 0.0) 
	{
		dptr->pvector = (atanf((vector1/vector2)) * (180.0/PI));
	}
	else if(vector2 < 0.0)
	{
		dptr->pvector = (atanf((vector1/vector2)) * (180.0/PI)) + 180.0;
	}
	else
	{
		dptr->pvector = (atanf((vector1/vector2)) * (180.0/PI))+ 360.0; 
	}	
	
	
	return;
}

void record_initialiser(data *dptr, params *pptr, int index)
{
	int cell;
	
	for(cell=0; cell<pptr->num_HD_cells; cell++)
	{
		dptr->initialiser_location_time[cell][index] = dptr->initialiser[cell];
	}
	
	return;
}


void fill_comb_buffer(data *dptr, params *pptr, cArray **bptr, int cell)
{
	bptr[cell]->array[bptr[cell]->index++] = dptr->rates_COMB[cell];
	
	if(bptr[cell]->index == pptr->conduction_buffer_size)
	{
		bptr[cell]->index = 0;
	}
	
	return;
}

inline float read_comb_buffer(cArray **bptr, int cell)
{
	return bptr[cell]->array[bptr[cell]->index];
}

void fill_hd_buffer(data *dptr, params *pptr, cArray **bxptr, int cell)
{
	bxptr[cell]->array[bxptr[cell]->index++] = dptr->rates_HD[cell];
	
	if(bxptr[cell]->index == pptr->conduction_buffer_size)
	{
		bxptr[cell]->index = 0;
	}
	
	return;
}

inline float read_hd_buffer(cArray **bxptr, int cell)
{
	return bxptr[cell]->array[bxptr[cell]->index];
}

