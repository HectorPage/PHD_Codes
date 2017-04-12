#ifndef __CALCULATE_INCLUDED__
#define __CALCULATE_INCLUDED__

#include "data.h"

extern void zero_states(data *dptr, params *pptr, cArray **bptr, connect *cptr);

extern void record_previous_states(data *dptr, params *pptr, connect *cptr);

extern void set_input(data *dptr, params *pptr);
extern void zero_input(data *dptr, params *pptr);


extern void calculate_activations_excitatory(data *dptr, params *pptr, connect *cptr, cArray **bptr, int timestep);
extern void calculate_rates_excitatory(data *dptr, params *pptr);

extern void record_activations_rates(data *dptr, params *pptr, int index);

extern void update_RC_weights(data *dptr, connect *cptr, params *pptr, cArray **bptr);

extern void record_input(data *dptr, params *pptr, int index);

extern void normalise_RC_weights(connect *cptr, params *pptr);

extern void fill_RC_buffer(data *dptr, params *pptr, cArray **bptr, int cell);
extern float read_RC_buffer(cArray **bptr, int cell);

extern void overlay_symmetrical(params *pptr, data *dptr, connect *cptr);

extern void calculate_weight_vector(data *dptr, params *pptr, connect *cptr, int index);

extern void log_weight_vector(data * dptr, params *pptr, connect *cptr, int index);


#endif