#ifndef __DATA_INCLUDED__
#define __DATA_INCLUDED__
#define PI 3.141592654


typedef struct {
	
	float time;
	
	int timesteps;
	
	int num_HD_cells;
	
	int num_COMB_cells;
	
	int num_ROT_cells;
	
	float phi_HDCOMB;
	
	float phi_COMBHD;
	
	float phi_ROT;
	
	float tau_HD;
	float tau_COMB;
	
	float alpha_HD;
	
	float beta_HD;
	
	float alpha_COMB;
	
	float beta_COMB;
	
	float timestep_size;
	
	float test_time;
	int test_timesteps;
	
	float pausetime;
	int pausetimesteps;
	
	float rotation_time;
	int rotation_timesteps;
	
	float NonOffTrainTime;
	int NonOffTrainTimesteps;
	
	float extern_inhib;
	float global_inhibition;
	float COMB_inhibition;
	
	float initialiser_str;
	float initialiser_loc;
	float initialiser_sigma;
		
	float conduction_delay;
	
	float velocity;
	
	int conduction_buffer_size;
	
	int numThreads;
	
	int introtest;
	
	float lrate;
	
	int training_loops;
		
	
}params;

typedef struct {
	
	float *activations_HD, *prev_activations_HD;

	float *rates_HD, *prev_rates_HD;
	
	float *activations_COMB, *prev_activations_COMB;
	
	float *rates_COMB, *prev_rates_COMB;

	float *rates_ROT, *prev_rates_ROT;
	
	float *initialiser;
		
	float **rates_HD_time;
		
	float **activations_HD_time;
	
	float **rates_COMB_time;
	
	float **activations_COMB_time;
	
	float **rates_ROT_time;
	
	float *favoured_view;
	
	float *time;
	
	float sumsq_HD_COMB;
	float sumsq_COMB_HD;
		
	float pvector;
	
	float **initialiser_location_time;
	
}data;

typedef struct {
	
	int num_HD_COMB_connections;
	int **HD_COMB_connections;
	//int HD_COMB_connections[1000][16]; 
	
	
	float **HD_COMB_weights;
	//float HD_COMB_weights[1000][16];
	float **prev_HD_COMB_weights;
	//float prev_HD_COMB_weights[1000][16];
	
	
	float **COMB_HD_weights;
	//float COMB_HD_weights[500][1000];
	float **prev_COMB_HD_weights;
	//float prev_COMB_HD_weights[500][1000];
	
	int inhib_random, excite_random;
	
	
	float **ROT_COMB_weights;
	//float ROT_COMB_weights[1000];
	float **prev_ROT_COMB_weights;
	//float prev_ROT_COMB_weights[1000];
	
	unsigned long int random_seed;	
	
	float **full_HDtoCOMB;
	
	
}connect;

typedef struct {

	int index;
	float *array;

}cArray;


#endif