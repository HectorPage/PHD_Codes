#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "calculate.h"
#include "data.h"

void zero_states(data *dptr, params *pptr)
{
	int cell;
	
	for(cell=0; cell<pptr->num_E1_cells; cell++)
	{
		dptr->activations_E1[cell] = 0.0;
		dptr->prev_activations_E1[cell] = 0.0;
		dptr->rates_E1[cell] = 0.0;
		dptr->prev_rates_E1[cell] = 0.0;
	}
	
		return;
}


void set_vis1(data *dptr, params *pptr)
{
    int cell;
    float distance, distance1, distance2;
    
    for(cell=0; cell<pptr->num_E1_cells; cell++)
    {
        distance1 = fabs(pptr->vis1_loc - dptr->favoured_view[cell]);
        distance2 = fabs(360.0 - distance1);
        
        if(distance1<=distance2)
            distance = distance1;
        else
            distance = distance2;
        
        dptr->visual_input1[cell] = pptr->vis1_str * exp(-0.5 * (distance/pptr->vis1_sigma) * (distance/pptr->vis1_sigma));
    }
}


void zero_vis1(data *dptr, params *pptr)
{
	int cell;
	
	for(cell=0; cell<pptr->num_E1_cells; cell++)
	{
		dptr->visual_input1[cell] = 0.0;
	}
	
	return;
}

void set_visring(data *dptr, params *pptr)
{
    int cell;
    float distance, distance1, distance2;
    
    for(cell=0; cell<pptr->num_E1_cells; cell++)
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
	
	for(cell=0; cell<pptr->num_E1_cells; cell++)
	{
		dptr->rates_visring[cell] = 0.0;
	}
	
	return;
}




void record_previous_states(data *dptr, params *pptr, connect *cptr)
{
	memcpy(dptr->prev_activations_E1, dptr->activations_E1, pptr->num_E1_cells*sizeof(float));
		
	memcpy(dptr->prev_rates_E1, dptr->rates_E1, pptr->num_E1_cells*sizeof(float));
	
	memcpy(dptr->prev_rates_visring, dptr->rates_visring, pptr->num_E1_cells*sizeof(float));
	
	memcpy(&(cptr->prev_weights_visring[0][0]), &(cptr->weights_visring[0][0]), pptr->num_E1_cells*pptr->num_E1_cells*sizeof(float));
	
	return;
}



void calculate_activations_excitatory(data *dptr, params *pptr, connect *cptr, int timestep)
{
	int post_cell, synapse, connection;
	float E1_coefficient;
	float inhibitory, recurrent;
	float RC1_scale;
	float visring;
	float scale_visring;
	
	E1_coefficient = pptr->timestep_size/pptr->tau_E1;
		
	RC1_scale = pptr->phi_RC1/(float)cptr->num_RC1_connections;
	
	scale_visring = pptr->phi_visring/(float)pptr->num_E1_cells;
	
	inhibitory = 0.0;
	
	for(synapse = 0; synapse<pptr->num_E1_cells; synapse++)
	{
		inhibitory += pptr->global_inhibition * dptr->prev_rates_E1[synapse];
	}
	
	for(post_cell=0; post_cell<pptr->num_E1_cells; post_cell++)
	{
		recurrent = 0.0;
		visring = 0.0;
		
		for(connection=0; connection<cptr->num_RC1_connections; connection++)
		{
			synapse = cptr->RC1_connections[post_cell][connection];
			
			recurrent += cptr->RC1_weights[post_cell][connection] * dptr->prev_rates_E1[synapse];
		}
		
		for(connection=0; connection<pptr->num_E1_cells; connection++)
		{
			visring += cptr->weights_visring[post_cell][connection] * dptr->prev_rates_visring[connection];
		}
		
		dptr->i_total[post_cell] = (dptr->sumsq * scale_visring * visring);
		
		dptr->activations_E1[post_cell] = ((1.0 - E1_coefficient) * dptr->prev_activations_E1[post_cell])
		+ (E1_coefficient * (dptr->sumsq * RC1_scale * recurrent))
		+ (E1_coefficient * dptr->visual_input1[post_cell])
		+ (E1_coefficient * dptr->sumsq * scale_visring * visring)
		- (E1_coefficient * inhibitory)
		- (E1_coefficient * pptr->extern_inhib);
	}
		
	return;
}

void calculate_pvector_vis1(data *dptr, params *pptr)
{
	int cell;
	float vector1 = 0.0;
	float vector2 = 0.0;
	
	for(cell=0; cell<pptr->num_E1_cells; cell++)
	{
		vector1+= dptr->rates_E1[cell] * dptr->favoured_view[cell];
		vector2+= dptr->rates_E1[cell];
	}
		dptr->pvector = vector1/vector2; 
		
	return;
}


void calculate_rates_excitatory(data *dptr, params *pptr)
{
	int cell;
	
	for(cell=0; cell<pptr->num_E1_cells; cell++)
	{
		dptr->rates_E1[cell] = 1.0 / (1.0 + exp(-2.0 * pptr->beta_E1 * (dptr->activations_E1[cell] - pptr->alpha_E1)));
	}
	
	return;
}

void record_activations_rates(data *dptr, params *pptr, int index)
{
	int cell;
	
	for(cell=0; cell<pptr->num_E1_cells; cell++)
	{
		dptr->activations_E1_time[cell][index] = dptr->activations_E1[cell];
		dptr->rates_E1_time[cell][index] = dptr->rates_E1[cell];
	}
	
	return;
}

void record_vis1(data *dptr, params *pptr, int index)
{
	int cell;
	
	for(cell=0; cell<pptr->num_E1_cells; cell++)
	{
		dptr->vis1_location_time[cell][index] = dptr->visual_input1[cell];
	}
	
	return;
}

void record_visring(data *dptr, params *pptr, int index)
{
	int cell;
		
	for(cell=0; cell<pptr->num_E1_cells; cell++)
	{
		dptr->visring_location_time[cell][index] = dptr->rates_visring[cell];
	}
		
	return;
}

void zero_averages(data *dptr, params *pptr)
{
	int cell;
	
	for(cell=0; cell<pptr->num_E1_cells; cell++)
	{
		dptr->E1_rates_average[cell] = 0.0;
		dptr->E1_rates_cumul[cell] = 0.0;
		
	}
	
	dptr->E1_rates_counter = 0;
	
	return;
}

void update_visring_weights(data *dptr, connect *cptr, params *pptr)
{
	
	int cell,synapse;
	
	
	for(cell=0;cell<pptr->num_E1_cells;cell++)
	{
		for(synapse=0;synapse<pptr->num_E1_cells;synapse++)
		{
			if(pptr->LTD)
			{
				cptr->weights_visring[cell][synapse] = cptr->prev_weights_visring[cell][synapse] + (pptr->lrate * pptr->timestep_size * (dptr->prev_rates_E1[cell] - dptr->E1_rates_average[cell]) * dptr->prev_rates_visring[synapse]);
				if(cptr->weights_visring[cell][synapse]<0.0)
				cptr->weights_visring[cell][synapse] = 0.0;
			}
			
			else
			{
				cptr->weights_visring[cell][synapse] = cptr->prev_weights_visring[cell][synapse] + (pptr->lrate * pptr->timestep_size * dptr->prev_rates_E1[cell] * dptr->prev_rates_visring[synapse]);
			}	
					
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
	
	for(post_cell = 0; post_cell<pptr->num_E1_cells; post_cell++)
	{
		sumsq = 0.0;
		
		for(pre_cell=0;pre_cell<pptr->num_E1_cells;pre_cell++)
		{
			sumsq+= cptr->weights_visring[post_cell][pre_cell]*cptr->weights_visring[post_cell][pre_cell];
		}
		
		if(sumsq==0.0)
			continue;
		
		sumsq = sqrt(sumsq);
		
		for(pre_cell=0;pre_cell<pptr->num_E1_cells;pre_cell++)
		{
			cptr->weights_visring[post_cell][pre_cell] = cptr->weights_visring[post_cell][pre_cell]/sumsq;
		}
	}
	
	return;
	
	}
			
	

void calculate_E1_rates_average(data *dptr, params *pptr)
{
	int cell;
	
	dptr->E1_rates_counter++;
	
	for(cell=0; cell<pptr->num_E1_cells; cell++)
	{
		dptr->E1_rates_cumul[cell] += dptr->rates_E1[cell];
		
		if(dptr->E1_rates_counter == 0)
			continue;
		
		dptr->E1_rates_average[cell] = dptr->E1_rates_cumul[cell]/(float)dptr->E1_rates_counter;
	}
	
	
	
	return;
}
	