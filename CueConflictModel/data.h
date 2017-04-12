#ifndef __DATA_INCLUDED__
#define __DATA_INCLUDED__


typedef struct {
	
	int timesteps;
	
	float time;
	
	int num_E1_cells;
		
	float phi_RC1;
	
	float tau_E1;
	
	float alpha_E1;
	
	float beta_E1;
	
	float timestep_size;
	
	float input_time;
	int input_timesteps;
	
	float extern_inhib;
	
	float global_inhibition;
	
	float vis1_str;
	float vis1_loc;
	float vis1_sigma;
	
	float visring_str;
	float visring_loc;
	float visring_sigma;
	
	float sigma_RC;
	float sigma_visring;
	
	float phi_visring;
	
	float lrate;
	
	int learning;
	int normalise;
	int LTD;
	
	int weight_dump_length;
	
	
}params;

typedef struct {
	
	float *activations_E1, *prev_activations_E1;
	
	float *rates_visring, *prev_rates_visring;

	float *rates_E1, *prev_rates_E1;
	
	float *visual_input1;
		
	float **rates_E1_time;
		
	float **activations_E1_time;
	
	float *favoured_view;
	
	float *time;
	
	float sumsq;
	
	float pvector;
	
	float **vis1_location_time;
	
	float **visring_location_time;
	
	float *E1_rates_average;
	
	float *E1_rates_cumul;
	
	int E1_rates_counter;
	
	float *i_total;
	
	
}data;

typedef struct {
	
	int num_RC1_connections;
	
	int **RC1_connections;
	float **RC1_weights;
	
	float **weights_visring;
	float **prev_weights_visring;
	
	int inhib_random, excite_random;
	
	unsigned long int random_seed;	
	
	
}connect;


#endif