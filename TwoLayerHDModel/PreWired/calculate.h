#ifndef __CALCULATE_INCLUDED__
#define __CALCULATE_INCLUDED__

#include "data.h"

extern void zero_states(data *dptr, params *pptr, cArray **bptr, cArray **bxptr);

extern void zero_COMB(data *dptr, params *pptr);

extern void record_previous_states(data *dptr, params *pptr, connect *cptr);

extern void set_initialiser(data *dptr, params *pptr);
extern void zero_initialiser(data *dptr, params *pptr);

extern void calculate_activations_excitatory(data *dptr, params *pptr, connect *cptr, cArray **bptr, cArray **bxptr, int timestep);
extern void calculate_rates_excitatory(data *dptr, params *pptr);

extern void record_activations_rates(data *dptr, params *pptr, int index);

extern void calculate_pvector_initialiser(data *dptr, params *pptr);

extern void record_initialiser(data *dptr, params *pptr, int index);

extern void fill_comb_buffer(data *dptr, params *pptr, cArray **bptr, int cell);

extern float read_comb_buffer(cArray **bptr, int cell);

extern void fill_hd_buffer(data *dptr, params *pptr, cArray **bxptr, int cell);

extern float read_hd_buffer(cArray **bxptr, int cell);


#endif