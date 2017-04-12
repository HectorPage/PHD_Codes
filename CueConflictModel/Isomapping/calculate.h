#ifndef __CALCULATE_INCLUDED__
#define __CALCULATE_INCLUDED__

#include "data.h"

extern void zero_states(data *dptr, params *pptr);

extern void record_previous_states(data *dptr, params *pptr, connect *cptr);

extern void set_pi(data *dptr, params *pptr);
extern void zero_pi(data *dptr, params *pptr);

extern void zero_visring(data *dptr, params *pptr);

extern void calculate_activations_excitatory(data *dptr, params *pptr, connect *cptr, int timestep);
extern void calculate_rates_excitatory(data *dptr, params *pptr);

extern void calculate_activations_inhibitory(data *dptr, params *pptr, connect *cptr);
extern void calculate_rates_inhibitory(data *dptr, params *pptr);

extern void record_activations_rates(data *dptr, params *pptr, int index);

extern void calculate_pvector(data *dptr, params *pptr);

extern void set_visring(data *dptr, params *pptr);

extern void update_visring_weights(data *dptr, connect *cptr, params *pptr);

extern void record_pi(data *dptr, params *pptr, int index);

extern void record_visring(data *dptr, params *pptr, int index);

extern void normalise_visring_weights(connect *cptr, params *pptr);

extern void calculate_weight_vector(data *dptr, params *pptr, connect *cptr, int index);

extern void initialise_pi(data *dptr, params *pptr);






#endif