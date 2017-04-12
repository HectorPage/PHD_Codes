#ifndef __DATA_INCLUDED__
#define __DATA_INCLUDED__
#define PI 3.141592654


typedef struct {
	
	int timesteps;
	
	float time;
	
	int num_HD_cells;
	
	float phi_RC1;
	
	float tau_HD;

	float alpha_HD;
	
	float beta_HD;
	
	float timestep_size;
	
	float input_time;
	int input_timesteps;
	float conflict_timesteps;
	int conflict_time;
	float pi_time;
	int pi_timesteps;
	
	float extern_inhib;
	float global_inhibition;
	
	float pi_str;
	float pi_loc;
	float pi_sigma;
	
	float visring_str;
	float visring_loc;
	float visring_sigma;
	
	float sigma_RC;
	float sigma_visring;
	
	float phi_visring;
	
	float lrate;
	int normalise;
	
	float velocity;
	
	int threads;

	
}params;

typedef struct {
	
	float *activations_HD, *prev_activations_HD;
	
	float *rates_visring, *prev_rates_visring;

	float *rates_HD, *prev_rates_HD;
	
	float *pi_input;
		
	float **rates_HD_time;
		
	float **activations_HD_time;	

	float *favoured_view;
	
	float *time;
	
	float sumsq_vis;
	float sumsq_RC;
	float pvector;
	
	float **pi_location_time;
	
	float **visring_location_time;
	
	float *HD_rates_average;
	
	float *HD_rates_cumul;
	
	int HD_rates_counter;
	
	float *weight_vector;
	
	float *HD_weight_vectors_init;
	
	float *HD_weight_vectors_intermediate;
	
	float *HD_weight_vectors_final;
	
	
	
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