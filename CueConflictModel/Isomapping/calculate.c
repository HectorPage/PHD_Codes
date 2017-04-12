#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "calculate.h"
#include "data.h"

void zero_states(data *dptr, params *pptr)
{
	int cell;

	for(cell=0; cell<pptr->num_HD_cells; cell++)
	{
		dptr->activations_HD[cell] = 0.0;
		dptr->prev_activations_HD[cell] = 0.0;
		dptr->rates_HD[cell] = 0.0;
		dptr->prev_rates_HD[cell] = 0.0;
	}
		return;
}

void initialise_pi(data *dptr, params *pptr)
{
    int cell;
    float distance, distancHD, distance2;
    
#pragma omp for private(cell, distancHD, distance2, distance)
    for(cell=0; cell<pptr->num_HD_cells; cell++)
    {
        distancHD = fabs(pptr->pi_loc - dptr->favoured_view[cell]);
        distance2 = fabs(360.0 - distancHD);
        
        if(distancHD<=distance2)
            distance = distancHD;
        else
            distance = distance2;
        
        dptr->pi_input[cell] = 10.0 * exp(-0.5 * (distance/pptr->pi_sigma) * (distance/pptr->pi_sigma));
    }
}




void set_pi(data *dptr, params *pptr)
{
    int cell;
    float distance, distancHD, distance2;
    
#pragma omp for private(cell, distancHD, distance2, distance)
    for(cell=0; cell<pptr->num_HD_cells; cell++)
    {
        distancHD = fabs(pptr->pi_loc - dptr->favoured_view[cell]);
        distance2 = fabs(360.0 - distancHD);
        
        if(distancHD<=distance2)
            distance = distancHD;
        else
            distance = distance2;
        
        dptr->pi_input[cell] = pptr->pi_str * exp(-0.5 * (distance/pptr->pi_sigma) * (distance/pptr->pi_sigma));
    }
}


void zero_pi(data *dptr, params *pptr)
{
	int cell;
#pragma omp for private(cell)
	for(cell=0; cell<pptr->num_HD_cells; cell++)
	{
		dptr->pi_input[cell] = 0.0;
	}
	
	return;
}

void set_visring(data *dptr, params *pptr)
{
    int cell;
    float distance, distance1, distance2;
    
#pragma omp for private(cell, distance, distance1, distance2) 
    for(cell=0; cell<pptr->num_HD_cells; cell++)
    {
        distance1 = fabs(pptr->visring_loc - dptr->favoured_view[cell]);
        distance2 = fabs(360.0 - distance1);
        
        if(distance1<=distance2)
            distance = distance1;
        else
            distance = distance2;
        
        dptr->rates_visring[cell] = pptr->visring_str * exp(-0.5 * (distance/pptr->visring_sigma) * (distance/pptr->visring_sigma));
    }
}



void zero_visring(data *dptr, params *pptr)
{
	int cell;
#pragma omp for private(cell)
	for(cell=0; cell<pptr->num_HD_cells; cell++)
	{
		dptr->rates_visring[cell] = 0.0;
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
			memcpy(&(dptr->prev_rates_visring[0]), &(dptr->rates_visring[0]), pptr->num_HD_cells*sizeof(float));
		}
#pragma omp section
		{
			memcpy(&(cptr->prev_weights_visring[0][0]), &(cptr->weights_visring[0][0]), pptr->num_HD_cells*pptr->num_HD_cells*sizeof(float));
		}
	}
	return;
}



void calculate_activations_excitatory(data *dptr, params *pptr, connect *cptr, int timestep)
{
	int post_cell, synapse, connection;
	float HD_coefficient;
	float HD_inhibitory, recurrent;
	float RC1_scale;
	float visring;
	float scale_visring;
	
	
	HD_coefficient = pptr->timestep_size/pptr->tau_HD;
	
		
	RC1_scale = pptr->phi_RC1/(float)cptr->num_RC1_connections;
	
	scale_visring = pptr->phi_visring/(float)pptr->num_HD_cells;
	
	/*Calculating Activation for HD cells*/
	
	HD_inhibitory = 0.0;
	
	for(synapse = 0; synapse<pptr->num_HD_cells; synapse++)
	{
		HD_inhibitory += pptr->global_inhibition * dptr->prev_rates_HD[synapse];
	}
#pragma omp for private(post_cell, visring, recurrent, connection)
	for(post_cell=0; post_cell<pptr->num_HD_cells; post_cell++)
	{
		recurrent = 0.0;
		visring = 0.0;
		
		for(connection=0; connection<cptr->num_RC1_connections; connection++)
		{
			synapse = cptr->RC1_connections[post_cell][connection];
			
			recurrent += cptr->RC1_weights[post_cell][connection] * dptr->prev_rates_HD[synapse];
		}
		
		for(connection=0; connection<pptr->num_HD_cells; connection++)
		{
			visring += cptr->weights_visring[post_cell][connection] * dptr->prev_rates_visring[connection];
		}
		
		dptr->activations_HD[post_cell] = ((1.0 - HD_coefficient) * dptr->prev_activations_HD[post_cell])
		+ (HD_coefficient * (dptr->sumsq_RC * RC1_scale * recurrent))
		+ (HD_coefficient * dptr->pi_input[post_cell])
		+ (HD_coefficient * scale_visring * visring)
		- (HD_coefficient * HD_inhibitory)
		- (HD_coefficient * pptr->extern_inhib);
		

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
		
	return;
}

void calculate_pvector(data *dptr, params *pptr)
{
	
	int cell;
	float vector1 = 0.0;
	float vector2 = 0.0;
	
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

void calculate_weight_vector(data *dptr, params *pptr, connect *cptr, int index)
{
	
	int cell;
	int connection;
	float vector1 = 0.0;
	float vector2 = 0.0;
	
	for(cell=0; cell<pptr->num_HD_cells; cell++)
	{
		vector1 = 0.0;
		vector2 = 0.0;
		
		for(connection=0; connection<pptr->num_HD_cells; connection++)
		{
			
		vector1+= cptr->weights_visring[cell][connection] * (sinf(dptr->favoured_view[connection]*(PI/180.0)));	/* Might be favoured view[cell] here???*/
		vector2+= cptr->weights_visring[cell][connection] * (cosf(dptr->favoured_view[connection]*(PI/180.0)));
			
		}
		
		if (vector1 > 0.0 && vector2 > 0.0) 
		{
			dptr->weight_vector[cell] = (atanf((vector1/vector2)) * (180.0/PI));
		}
		else if(vector2 < 0.0)
		{
			dptr->weight_vector[cell] = (atanf((vector1/vector2)) * (180.0/PI)) + 180.0;
		}
		else
		{
			dptr->weight_vector[cell] = (atanf((vector1/vector2)) * (180.0/PI)) + 360.0;
		}	
		
		switch (index) {
			case 0:
				dptr->HD_weight_vectors_init[cell] = dptr->weight_vector[cell];
				break;
			case 1:
				dptr->HD_weight_vectors_intermediate[cell] = dptr->weight_vector[cell];
				break;
			case 2:
				dptr->HD_weight_vectors_final[cell] = dptr->weight_vector[cell];
		}
			
	}
	return;
}

void record_pi(data *dptr, params *pptr, int index)
{
	int cell;
	
	for(cell=0; cell<pptr->num_HD_cells; cell++)
	{
		dptr->pi_location_time[cell][index] = dptr->pi_input[cell];
	}
	
	return;
}

void record_visring(data *dptr, params *pptr, int index)
{
	int cell;
		
	for(cell=0; cell<pptr->num_HD_cells; cell++)
	{
		dptr->visring_location_time[cell][index] = dptr->rates_visring[cell];
	}
		
	return;
}

void update_visring_weights(data *dptr, connect *cptr, params *pptr)
{
	
	int cell,synapse;
	
#pragma omp for private(cell, synapse)	
	for(cell=0;cell<pptr->num_HD_cells;cell++)
	{
		for(synapse=0;synapse<pptr->num_HD_cells;synapse++)
		{
				cptr->weights_visring[cell][synapse] = cptr->prev_weights_visring[cell][synapse] + (pptr->lrate * pptr->timestep_size * dptr->prev_rates_HD[cell] * dptr->prev_rates_visring[synapse]);
						
		}

	}

	
	if(pptr->normalise)
	{
		normalise_visring_weights(cptr, pptr);
	}
	
	return;
}

void normalise_visring_weights(connect *cptr, params *pptr)
{
	int post_cell, pre_cell;
	
	float sumsq;
#pragma omp for private(post_cell, pre_cell, sumsq)	
	for(post_cell = 0; post_cell<pptr->num_HD_cells; post_cell++)
	{
		sumsq = 0.0;
		
		for(pre_cell=0;pre_cell<pptr->num_HD_cells;pre_cell++)
		{
			sumsq+= cptr->weights_visring[post_cell][pre_cell]*cptr->weights_visring[post_cell][pre_cell];
		}
		
		if(sumsq==0.0)
			continue;
		
		sumsq = sqrt(sumsq);
		
		for(pre_cell=0;pre_cell<pptr->num_HD_cells;pre_cell++)
		{
			cptr->weights_visring[post_cell][pre_cell] = cptr->weights_visring[post_cell][pre_cell]/sumsq;
		}
	}
	
	return;
	
	}
		