#ifndef __CALCULATE_INCLUDED__
#define __CALCULATE_INCLUDED__

#include "data.h"

#define eps 1e-7
extern void zero_states(data *dptr, params *pptr, cArray **combptr, cArray **hdptr);

extern void zero_COMB(data *dptr, params *pptr);

extern void record_previous_states(data *dptr, params *pptr, connect *cptr, cArray **hdptr);

extern void set_initialiser(data *dptr, params *pptr);
extern void zero_initialiser(data *dptr, params *pptr);

extern void calculate_activations_excitatory(data *dptr, params *pptr, connect *cptr, cArray **combptr, cArray **hdptr);
extern void calculate_rates_excitatory(data *dptr, params *pptr);

extern void record_activations_rates(data *dptr, params *pptr, int index);

extern void record_initialiser(data *dptr, params *pptr, int index);

extern void fill_comb_buffer(data *dptr, params *pptr, cArray **combptr, int cell);

extern float read_comb_buffer(cArray **combptr, int cell);

extern void fill_hd_buffer(data *dptr, params *pptr, cArray **hdptr, int cell);

extern float read_hd_buffer(cArray **hdptr, int cell);

extern void update_HDtoCOMB_weights(data *dptr, connect *cptr, params *pptr, cArray **hdptr);

extern void update_COMBtoHD_weights(data *dptr, connect *cptr, params *pptr, cArray **combptr);

extern void normalise_HDtoCOMB_weights(connect *cptr, params *pptr);

extern void normalise_COMBtoHD_weights(connect *cptr, params *pptr);

extern void normalise_ROTtoCOMB_weights(connect *cptr, params *pptr);

extern void update_ROTtoCOMB_weights(data *dptr, connect *cptr, params *pptr);

extern void zero_previous_states(params *pptr, connect *cptr, data *dptr);

extern void calculate_pvector(data *dptr, params *pptr);

#endif