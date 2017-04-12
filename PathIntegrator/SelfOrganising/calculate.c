#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "calculate.h"
#include "data.h"

void zero_states(data *dptr, params *pptr, cArray **bptr, connect *cptr)
{
	int cell;
	int connection;
	
	for(cell=0; cell<pptr->num_E1_cells; cell++)
	{
		dptr->activations_E1[cell] = 0.0;
		dptr->prev_activations_E1[cell] = 0.0;
		dptr->rates_E1[cell] = 0.0;
		dptr->prev_rates_E1[cell] = 0.0;
	}
	
	for(cell=0; cell<pptr->conduction_buffer_size; cell++)
	{
		for(connection=0; connection<pptr->num_E1_cells; connection++) //Is this num RC connections?
		{
			bptr[connection]->array[cell] = 0.0;			
		}
	}
	
	
		return;
}


void set_input(data *dptr, params *pptr)
{
    int cell;
    float distance, distance1, distance2;
    
#pragma omp for private(cell, distance, distance1, distance2)
    for(cell=0; cell<pptr->num_E1_cells; cell++)
    {
        distance1 = fabs(pptr->input_loc - dptr->favoured_view[cell]);
        distance2 = fabs(360.0 - distance1);
        
        if(distance1<=distance2)
            distance = distance1;
        else
            distance = distance2;
        
        dptr->input[cell] = pptr->input_str * exp(-0.5 * (distance/pptr->input_sigma) * (distance/pptr->input_sigma));
    }
}


void zero_input(data *dptr, params *pptr)
{
	int cell;
#pragma omp for private(cell)
	for(cell=0; cell<pptr->num_E1_cells; cell++)
	{
		dptr->input[cell] = 0.0;
	}
	
	return;
}



void record_previous_states(data *dptr, params *pptr, connect *cptr)
{
#pragma omp sections
	{
#pragma omp section
		{
			memcpy(dptr->prev_activations_E1, dptr->activations_E1, pptr->num_E1_cells*sizeof(float));
		}
#pragma omp section
		{
			memcpy(dptr->prev_rates_E1, dptr->rates_E1, pptr->num_E1_cells*sizeof(float));
		}
#pragma omp section
		{
			memcpy(&(cptr->prev_weights_RC1[0][0]), &(cptr->RC1_weights[0][0]), pptr->num_E1_cells*cptr->num_RC1_connections*sizeof(float));
		}
	}
	return;
}



void calculate_activations_excitatory(data *dptr, params *pptr, connect *cptr, cArray **bptr, int timestep)
{
	int post_cell, synapse, connection;
	float E1_coefficient;
	float inhibitory, recurrent;
	float RC1_scale;

	
	int cArray_index;
	
	cArray_index = 0;
	
	E1_coefficient = pptr->timestep_size/pptr->tau_E1;
		
	if(pptr->testphase)
	{
		RC1_scale = pptr->phi_RC1_test/(float)cptr->num_RC1_connections;

	}
	else {
		RC1_scale = pptr->phi_RC1/(float)cptr->num_RC1_connections;
	}
	
	inhibitory = 0.0;
	
	for(synapse = 0; synapse<pptr->num_E1_cells; synapse++)
	{
		inhibitory += pptr->global_inhibition * dptr->prev_rates_E1[synapse];
	}
	

	
#pragma omp for private(post_cell, recurrent, connection, synapse)
	for(post_cell=0; post_cell<pptr->num_E1_cells; post_cell++)
	{
		recurrent = 0.0;
				
		for(connection=0; connection<cptr->num_RC1_connections; connection++)
		{
			synapse = cptr->RC1_connections[post_cell][connection];
			
			recurrent += cptr->RC1_weights[post_cell][synapse] * read_RC_buffer(bptr, synapse); //Could be RC1_weights[post_cell][connection]?
		}
		
		dptr->activations_E1[post_cell] = ((1.0 - E1_coefficient) * dptr->prev_activations_E1[post_cell])
		+ (E1_coefficient * (RC1_scale * recurrent))
		+ (E1_coefficient * dptr->input[post_cell])
		- (E1_coefficient * inhibitory)
		- (E1_coefficient * pptr->extern_inhib);
	}
		
	return;
}


void calculate_rates_excitatory(data *dptr, params *pptr)
{
	int cell;
	
	if(pptr->sigmoid)
	   {

#pragma omp for private(cell)
		   for(cell=0; cell<pptr->num_E1_cells; cell++)
			{
		
			dptr->rates_E1[cell] = 1.0 / (1.0 + exp(-2.0 * pptr->beta_E1 * (dptr->activations_E1[cell] - pptr->alpha_E1)));
	
			}
	   }
	   
	   else
	   {
#pragma omp for private(cell)
	   for(cell=0; cell<pptr->num_E1_cells; cell++)
	   {
		   
		   dptr->rates_E1[cell] = (float)tanh(dptr->activations_E1[cell]);
		   
		   if(dptr->rates_E1[cell] < 0.0)
			   dptr->rates_E1[cell] = 0.0;
		   
	   }
		   
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

void record_input(data *dptr, params *pptr, int index)
{
	int cell;
	
	for(cell=0; cell<pptr->num_E1_cells; cell++)
	{
		dptr->input_location_time[cell][index] = dptr->input[cell];
	}
	
	return;
}




void update_RC_weights(data *dptr, connect *cptr, params *pptr, cArray **bptr)
{
	
	int cell,synapse, presynaptic;
	
#pragma omp for private(cell, synapse, presynaptic)
	for(cell=0;cell<pptr->num_E1_cells;cell++)
	{
		for(synapse=0;synapse<cptr->num_RC1_connections;synapse++)
		{
			
			presynaptic = cptr->RC1_connections[cell][synapse];
			
			cptr->RC1_weights[cell][presynaptic] = cptr->prev_weights_RC1[cell][presynaptic] + (pptr->lrate * pptr->timestep_size * dptr->prev_rates_E1[cell] * read_RC_buffer(bptr, presynaptic));
			
			//cptr->RC1_weights[cell][synapse] = cptr->prev_weights_RC1[cell][synapse] + (pptr->lrate * pptr->timestep_size * dptr->prev_rates_E1[cell] * read_RC_buffer(bptr, presynaptic));
									
			/*if(presynaptic == cell)
			{
				cptr->RC1_weights[cell][synapse] = 0.0;
			}*/
			
			
			
		}

	}

	
	if(pptr->normalise)
	{
		normalise_RC_weights(cptr, pptr);
	}
	
	return;
}

void normalise_RC_weights(connect *cptr, params *pptr)
{
	int post_cell, synapse, presynaptic;
	
	float sumsq;
#pragma omp for private(post_cell, synapse, sumsq, presynaptic)	
	for(post_cell = 0; post_cell<pptr->num_E1_cells; post_cell++)
	{
		sumsq = 0.0;
		
		for(synapse=0;synapse<cptr->num_RC1_connections;synapse++)
		{
			presynaptic = cptr->RC1_connections[post_cell][synapse];
			
			sumsq+= cptr->RC1_weights[post_cell][presynaptic]*cptr->RC1_weights[post_cell][presynaptic];
		}
		
		if(sumsq==0.0)
			continue;
		
		sumsq = sqrt(sumsq);
		
		for(synapse=0;synapse<cptr->num_RC1_connections;synapse++)
		{
			presynaptic = cptr->RC1_connections[post_cell][synapse];

			cptr->RC1_weights[post_cell][presynaptic] = cptr->RC1_weights[post_cell][presynaptic]/sumsq;
		}
	}
	
	return;
	
}

void fill_RC_buffer(data *dptr, params *pptr, cArray **bptr, int cell)
{
	bptr[cell]->array[bptr[cell]->index++] = dptr->rates_E1[cell];
	
	if(bptr[cell]->index == pptr->conduction_buffer_size)
	{
		bptr[cell]->index = 0;
	}
	
	return;
	
}

inline float read_RC_buffer(cArray **bptr, int cell)
{
	return bptr[cell]->array[bptr[cell]->index];
}



void overlay_symmetrical(params *pptr, data *dptr, connect *cptr)
{
	int cell, connection, presynaptic;
	float increment, distance1, distance2, distance, extra;
	
	increment = 360.0/(float)cptr->num_RC1_connections;
	
#pragma omp for private(cell, connection, distance1, distance2, distance, presynaptic, extra)
	for(cell=0; cell<pptr->num_E1_cells; cell++)
	{
		
		for(connection=0; connection<cptr->num_RC1_connections; connection++)
		{

			presynaptic = cptr->RC1_connections[cell][connection];
			
			distance1 = fabs(dptr->favoured_view[presynaptic] - dptr->favoured_view[cell]);			
			distance2 = fabs(360.0 - distance1);
			
			if(distance1<=distance2)
				distance = distance1;
			else
			{
				distance = distance2;
			}
			
			extra = pptr->sym_str * (exp(-0.5 * (distance/pptr->sym_sigma) * (distance/pptr->sym_sigma)));
			
			cptr->RC1_weights[cell][presynaptic] += extra;	
			
		}
	}
	
	normalise_RC_weights(cptr, pptr);
	
#pragma omp single
	{
		printf("\n\nSymmetrical weights overlaid\n\n");
		fflush(stdout);
	}
	
	return;
	
}


void calculate_weight_vector(data *dptr, params *pptr, connect *cptr, int index)
{
	//This might not work for incomplete connectivity.....
	
	int postsynaptic;
	int presynaptic;
	float vector1 = 0.0;
	float vector2 = 0.0;
	
	for(presynaptic=0; presynaptic<cptr->num_RC1_connections; presynaptic++)
	{
		vector1 = 0.0;
		vector2 = 0.0;
		
		for(postsynaptic=0; postsynaptic<pptr->num_E1_cells; postsynaptic++)
		{

			vector1+= cptr->full_RC[postsynaptic][presynaptic] * (sinf(dptr->favoured_view[postsynaptic]*(PI/180.0)));
			vector2+= cptr->full_RC[postsynaptic][presynaptic] * (cosf(dptr->favoured_view[postsynaptic]*(PI/180.0)));
			
		}
		
		if (vector1 > 0.0 && vector2 > 0.0) 
		{
			dptr->weight_vector[presynaptic] = (atanf((vector1/vector2)) * (180.0/PI));
		}
		else if(vector2 < 0.0)
		{
			dptr->weight_vector[presynaptic] = (atanf((vector1/vector2)) * (180.0/PI)) + 180.0;
		}
		else
		{
			dptr->weight_vector[presynaptic] = (atanf((vector1/vector2)) * (180.0/PI)) + 360.0;
		}	
		
		switch (index) {
			case 1:
				dptr->weight_vector_intermediate[presynaptic] = dptr->weight_vector[presynaptic];
				break;
			case 2:
				dptr->weight_vector_final[presynaptic] = dptr->weight_vector[presynaptic];
		}
		
	}
	return;
}

void log_weight_vector(data * dptr, params *pptr, connect *cptr, int index)
{
	
	int postsynaptic;
	int presynaptic;
	float vector1 = 0.0;
	float vector2 = 0.0;
	
	for(presynaptic=0; presynaptic<cptr->num_RC1_connections; presynaptic++)
	{
		vector1 = 0.0;
		vector2 = 0.0;
		
		for(postsynaptic=0; postsynaptic<pptr->num_E1_cells; postsynaptic++)
		{
			
			vector1+= cptr->full_RC[postsynaptic][presynaptic] * (sinf(dptr->favoured_view[postsynaptic]*(PI/180.0)));
			vector2+= cptr->full_RC[postsynaptic][presynaptic] * (cosf(dptr->favoured_view[postsynaptic]*(PI/180.0)));
			
		}
		
		if (vector1 > 0.0 && vector2 > 0.0) 
		{
			dptr->wVec[presynaptic] = (atanf((vector1/vector2)) * (180.0/PI));
		}
		else if(vector2 < 0.0)
		{
			dptr->wVec[presynaptic] = (atanf((vector1/vector2)) * (180.0/PI)) + 180.0;
		}
		else
		{
			dptr->wVec[presynaptic] = (atanf((vector1/vector2)) * (180.0/PI)) + 360.0;
		}	
		
		dptr->wVec_time[presynaptic][index] = dptr->wVec[presynaptic];

	}
			
	
	return;
}

	

	