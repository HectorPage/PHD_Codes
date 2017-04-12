#ifndef __DATA_INCLUDED__
#define __DATA_INCLUDED__
#define PI 3.141592654


typedef struct {
	
	int timesteps;
	
	float time;
	
	int num_E1_cells;
		
	float phi_RC1;
	float phi_RC1_test;
	
	float tau_E1;
	
	float alpha_E1;
	
	float beta_E1;
	
	float timestep_size;
	
	float training_time;
	int training_timesteps;
	
	float extern_inhib;
	
	float global_inhibition;
	
	float input_str;
	float input_loc;
	float input_sigma;

	int normalise;
	
	float velocity;
	
	float conduction_delay;
	int conduction_buffer_size;
	
	int threads;
	
	int testphase;
	
	int sigmoid;
	
	int symtest;
	
	float sym_sigma;
	
	float sym_str;
	
}params;

typedef struct {
	
	float *activations_E1, *prev_activations_E1;

	float *rates_E1, *prev_rates_E1;
	
	float *input;
		
	float **rates_E1_time;
		
	float **activations_E1_time;
	
	float *favoured_view;
	
	float *time;
	
	float sumsq;
	
	float **input_location_time;

	float *E1_rates_average;
	
	float *E1_rates_cumul;
	
	int E1_rates_counter;
	
	float *weight_vector;
	
	float *weight_vector_intermediate;
	
	float *weight_vector_final;
	
	
}data;

typedef struct {
	
	int num_RC1_connections;
	
	int **RC1_connections;
	
	float **RC1_weights;
	float **prev_weights_RC1;
	
	float **full_RC;
		
	
	int inhib_random, excite_random;
	
	unsigned long int random_seed;	
	
	
}connect;

typedef struct {
	
	int index;
	float *array;
}cArray;


#endif