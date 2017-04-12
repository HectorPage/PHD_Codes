#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "calculate.h"
#include "data.h"

void zero_states(data *dptr, params *pptr, cArray **combptr, cArray **hdptr)
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
	
	
	for(connection=0; connection<pptr->num_COMB_cells; connection++)
	{
		for(cell=0; cell<pptr->conduction_buffer_size; cell++)
		{
			
			combptr[connection]->array[cell] = 0.0;
			
		}
	}
	
	printf("COMB buffer cleared\n");
	fflush(stdout);

	for(connection=0; connection<pptr->num_HD_cells; connection++)
	{
		
		for(cell=0; cell<pptr->conduction_buffer_size; cell++)
		{
			
			hdptr[connection]->array[cell] = 0.0;
			
		}
	}
	
	
	printf("HD buffer cleared\n");
	fflush(stdout);
	
		
	for(cell=0;cell<pptr->num_ROT_cells;cell++)
	{
		dptr->rates_ROT[cell] = 0.0;
	}
	
	printf("ROT rates zeroed\n");
	fflush(stdout);


	
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


void record_previous_states(data *dptr, params *pptr, connect *cptr, cArray **hdptr)
{
	
#pragma omp parallel sections
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
			
#pragma omp section
		{
			memcpy(&(cptr->prev_HD_COMB_weights[0][0]), &(cptr->HD_COMB_weights[0][0]), pptr->num_COMB_cells*cptr->num_HD_COMB_connections*sizeof(float));

			
		}
	
#pragma omp section
		{
			memcpy(&(cptr->prev_COMB_HD_weights[0][0]), &(cptr->COMB_HD_weights[0][0]), pptr->num_HD_cells*pptr->num_COMB_cells*sizeof(float));

			

		}

#pragma omp section
		{
			memcpy(&(dptr->prev_rates_ROT[0]), &(dptr->rates_ROT[0]), pptr->num_ROT_cells*sizeof(float));

			
		}	
#pragma omp section
		{
			memcpy(&(cptr->prev_ROT_COMB_weights[0][0]), &(cptr->ROT_COMB_weights[0][0]), pptr->num_COMB_cells*pptr->num_ROT_cells*sizeof(float));


		}
		
		
	}
	 
	 
}	



void calculate_activations_excitatory(data *dptr, params *pptr, connect *cptr, cArray **combptr, cArray **hdptr)
{
	int post_cell, synapse, connection, presynaptic_cell, ROTcount;
	float HD_coefficient, COMB_coefficient;
	float HD_inhibitory, COMB_inhibitory, ROT_input, COMB_input, HD_input_COMB;
	float HD_COMB_scale, COMB_HD_scale, ROT_COMB_scale;
	
	HD_coefficient = pptr->timestep_size/pptr->tau_HD;
	COMB_coefficient = pptr->timestep_size/pptr->tau_COMB;
		
	HD_COMB_scale = pptr->phi_HDCOMB/(float)cptr->num_HD_COMB_connections;
	COMB_HD_scale = pptr->phi_COMBHD/(float)pptr->num_COMB_cells;
	ROT_COMB_scale = pptr->phi_ROT/(float)pptr->num_ROT_cells; 
		
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
		- (HD_coefficient * HD_inhibitory *(1.0/pptr->num_HD_cells))
		- (HD_coefficient * pptr->extern_inhib);
		

	}
	
	/*Calculating Activation for COMB cells*/
	
		
	COMB_inhibitory = 0.0;
	
	for(synapse=0; synapse<pptr->num_COMB_cells; synapse++)
	{ 
		COMB_inhibitory += pptr->COMB_inhibition * dptr->prev_rates_COMB[synapse];
	}
	
#pragma omp for private(post_cell, HD_input_COMB, connection, presynaptic_cell, ROT_input, ROTcount)
	for(post_cell=0; post_cell<pptr->num_COMB_cells; post_cell++)
	{
		HD_input_COMB= 0.0;
		
		for(connection=0; connection<cptr->num_HD_COMB_connections; connection++)
		{
			presynaptic_cell = cptr->HD_COMB_connections[post_cell][connection];
			HD_input_COMB += cptr->HD_COMB_weights[post_cell][connection] * read_hd_buffer(hdptr, presynaptic_cell);
		}
			
		ROT_input = 0.0;
				
		for(ROTcount=0; ROTcount<pptr->num_ROT_cells; ROTcount++)
		{
			ROT_input += cptr->ROT_COMB_weights[post_cell][ROTcount] * dptr->prev_rates_ROT[ROTcount];
		}
		
		
	
		dptr->activations_COMB[post_cell] = ((1.0 - COMB_coefficient) * dptr->prev_activations_COMB[post_cell])
		+ ((HD_input_COMB * HD_COMB_scale) * COMB_coefficient)								
		+ ((ROT_input * ROT_COMB_scale) * COMB_coefficient)
		- (COMB_inhibitory * COMB_coefficient * (1.0/pptr->num_COMB_cells));
		
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
	
	for(cell=0; cell<pptr->num_ROT_cells; cell++)
	{
		dptr->rates_ROT_time[cell][index] = dptr->rates_ROT[cell];
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


void fill_comb_buffer(data *dptr, params *pptr, cArray **combptr, int cell)
{
	combptr[cell]->array[combptr[cell]->index++] = dptr->rates_COMB[cell];
	
	if(combptr[cell]->index == pptr->conduction_buffer_size)
	{
		combptr[cell]->index = 0;
	}
	
	return;
}

inline float read_comb_buffer(cArray **combptr, int cell)
{
	return combptr[cell]->array[combptr[cell]->index];
}

void fill_hd_buffer(data *dptr, params *pptr, cArray **hdptr, int cell)
{
	hdptr[cell]->array[hdptr[cell]->index++] = dptr->rates_HD[cell];
	
	if(hdptr[cell]->index == pptr->conduction_buffer_size)
	{
		hdptr[cell]->index = 0;
	}
	
	return;
}

inline float read_hd_buffer(cArray **hdptr, int cell)
{
	return hdptr[cell]->array[hdptr[cell]->index];
}

void update_HDtoCOMB_weights(data *dptr, connect *cptr, params *pptr, cArray **hdptr)
{
	int cell, synapse, presynaptic_cell;

	
#pragma omp for private(cell, synapse, presynaptic_cell)
	for(cell=0;cell<pptr->num_COMB_cells; cell++)
	{
		for(synapse=0;synapse<cptr->num_HD_COMB_connections;synapse++)
		{
			presynaptic_cell = cptr->HD_COMB_connections[cell][synapse];
			
			cptr->HD_COMB_weights[cell][synapse] = cptr->prev_HD_COMB_weights[cell][synapse] + (pptr->lrate * pptr->timestep_size * dptr->prev_rates_COMB[cell] * read_hd_buffer(hdptr, presynaptic_cell));
		}
	}
	
			
	normalise_HDtoCOMB_weights(cptr, pptr);
		
	return;
}
	
void update_COMBtoHD_weights(data *dptr, connect *cptr, params *pptr, cArray **combptr)
{
		int cell, synapse;
			
#pragma omp for private(cell, synapse)
		for(cell=0;cell<pptr->num_HD_cells; cell++)
		{
			for(synapse=0;synapse<pptr->num_COMB_cells;synapse++)
			{
				cptr->COMB_HD_weights[cell][synapse] = cptr->prev_COMB_HD_weights[cell][synapse] + (pptr->lrate * pptr->timestep_size * dptr->prev_rates_HD[cell] * read_comb_buffer(combptr, synapse));
			}
		}
	
	normalise_COMBtoHD_weights(cptr, pptr);
	
	
	return;
	
}

void normalise_HDtoCOMB_weights(connect *cptr, params *pptr)
{
	int post_cell, connection;
	
	float sumsq;
#pragma omp for private(post_cell,sumsq, connection)	
	for(post_cell = 0; post_cell<pptr->num_COMB_cells; post_cell++)
	{
		sumsq = 0.0;
		
		for(connection=0;connection<cptr->num_HD_COMB_connections;connection++)
		{
			sumsq+= cptr->HD_COMB_weights[post_cell][connection]*cptr->HD_COMB_weights[post_cell][connection];
		}
		
		if(fabs(sumsq)<eps)
			continue;
		
		sumsq = sqrt(sumsq);
		
		for(connection=0;connection<cptr->num_HD_COMB_connections;connection++)
		{			
			cptr->HD_COMB_weights[post_cell][connection] = cptr->HD_COMB_weights[post_cell][connection]/sumsq;
		}
	}
	
	return;
}

void normalise_COMBtoHD_weights(connect *cptr, params *pptr)
{
	int post_cell, pre_cell;
	
	float sumsq;
#pragma omp for private(post_cell, pre_cell, sumsq)	
	for(post_cell = 0; post_cell<pptr->num_HD_cells; post_cell++)
	{
		sumsq = 0.0;
		
		for(pre_cell=0;pre_cell<pptr->num_COMB_cells;pre_cell++)
		{
			sumsq+= cptr->COMB_HD_weights[post_cell][pre_cell]*cptr->COMB_HD_weights[post_cell][pre_cell];
		}
		
		if(fabs(sumsq)<eps)
			continue;
		
		sumsq = sqrt(sumsq);
		
		for(pre_cell=0;pre_cell<pptr->num_COMB_cells;pre_cell++)
		{
			cptr->COMB_HD_weights[post_cell][pre_cell] = cptr->COMB_HD_weights[post_cell][pre_cell]/sumsq;
		}
	}
	
	return;
}

void update_ROTtoCOMB_weights(data *dptr, connect *cptr, params *pptr)
{
	int cell;
	int connection;
	
#pragma omp for private(cell)
	for(cell=0;cell<pptr->num_COMB_cells; cell++)
	{
		for(connection=0;connection<pptr->num_ROT_cells; connection++)
		{
		
			cptr->ROT_COMB_weights[cell][connection] = cptr->prev_ROT_COMB_weights[cell][connection] + (pptr->lrate * pptr->timestep_size * dptr->prev_rates_COMB[cell] * dptr->prev_rates_ROT[connection]);
		}
	}
	
	
	normalise_ROTtoCOMB_weights(cptr, pptr);
	
	return;
}

void normalise_ROTtoCOMB_weights(connect *cptr, params *pptr)
{
	int post_cell;
	int pre_cell;
	
	float sumsq;

	for(post_cell = 0; post_cell<pptr->num_COMB_cells; post_cell++)
	{
		sumsq  = 0.0;
		
		for(pre_cell=0; pre_cell<pptr->num_ROT_cells; pre_cell++)
		{
			sumsq += cptr->ROT_COMB_weights[post_cell][pre_cell]*cptr->ROT_COMB_weights[post_cell][pre_cell];
		}
		
		if(fabs(sumsq)<eps)
			continue;
		
		sumsq = sqrt(sumsq);
		
		for(pre_cell=0; pre_cell<pptr->num_ROT_cells; pre_cell++)
		{
			cptr->ROT_COMB_weights[post_cell][pre_cell] = cptr->ROT_COMB_weights[post_cell][pre_cell]/sumsq;
		}

	}
		return;
}


void zero_previous_states(params *pptr, connect *cptr, data *dptr)
{
	int cell, connection;
		
	for(connection=0;connection<pptr->num_ROT_cells;connection++)
	{
		dptr->rates_ROT[connection] = 0.0;
	}
	
	for(cell=0;cell<pptr->num_HD_cells;cell++)
	{
		dptr->prev_activations_HD[cell] = 0.0;
		dptr->prev_rates_HD[cell] = 0.0;
		
		for(connection=0;connection<pptr->num_COMB_cells;connection++)
		{
			cptr->prev_COMB_HD_weights[cell][connection] = 0.0;
		}
	}
	
	for(cell=0;cell<pptr->num_COMB_cells;cell++)
	{
		dptr->prev_activations_COMB[cell] = 0.0;
		dptr->prev_rates_COMB[cell] = 0.0;
		
		for(connection=0; connection<pptr->num_ROT_cells; connection++)
		{
			cptr->prev_ROT_COMB_weights[cell][connection] = 0.0;
		}
				
		for(connection=0;connection<cptr->num_HD_COMB_connections;connection++)
		{
			cptr->prev_HD_COMB_weights[cell][connection] = 0.0;
		}
	}
	
	
	
	
	return;
}

void calculate_pvector(data *dptr, params *pptr)
{
	
	int cell;
	float vector1 = 0.0;
	float vector2 = 0.0;

#pragma omp for private(cell)
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

