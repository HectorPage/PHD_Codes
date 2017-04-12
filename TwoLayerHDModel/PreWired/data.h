#ifndef __DATA_INCLUDED__
#define __DATA_INCLUDED__
#define PI 3.141592654


typedef struct {
	
	int timesteps;
	
	float time;
	
	int num_HD_cells;
	
	int num_COMB_cells;
	
	int num_ROT_cells;
	
	int num_NOROT_cells;
	
	float phi_HDCOMB;
	
	float phi_COMBHD;
	
	float phi_ROT;
	
	float phi_NOROT;
	
	float tau_HD;
	float tau_COMB;
	
	float alpha_HD;
	
	float beta_HD;
	
	float alpha_COMB;
	
	float beta_COMB;
	
	float timestep_size;
	
	float input_time;
	int input_timesteps;
	
	float pausetime;
	int pausetimesteps;
	
	float rotation_time;
	int rotation_timesteps;
	
	float extern_inhib;
	float global_inhibition;
	float COMB_inhibition;
	
	float initialiser_str;
	float initialiser_loc;
	float initialiser_sigma;
	
	float sigma_COMB_HD;
	float sigma_HD_COMB;	
	
	float conduction_delay;
	
	float velocity;
	
	int conduction_buffer_size;
	
	int numThreads;
		
	
}params;

typedef struct {
	
	float *activations_HD, *prev_activations_HD;

	float *rates_HD, *prev_rates_HD;
	
	float *activations_COMB, *prev_activations_COMB;
	
	float *rates_COMB, *prev_rates_COMB;

	float rates_ROT;
	
	float rates_NOROT;
	
	float *initialiser;
		
	float **rates_HD_time;
		
	float **activations_HD_time;
	
	float **rates_COMB_time;
	
	float **activations_COMB_time;
	
	float *rates_ROT_time;
	
	float *rates_NOROT_time;
	
	float *favoured_view;
	
	float *time;
	
	float sumsq_HD_COMB;
	float sumsq_COMB_HD;
		
	float pvector;
	
	float **initialiser_location_time;
	
}data;

typedef struct {
	
	int **HD_COMB_connections;
	float **HD_COMB_weights;
	

	int **COMB_HD_connections;
	float **COMB_HD_weights;
	
	int inhib_random, excite_random;
	
	unsigned long int random_seed;	
	
	
}connect;

typedef struct {

	int index;
	float *array;

}cArray;


#endif