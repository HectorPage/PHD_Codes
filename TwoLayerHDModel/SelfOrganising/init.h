#ifndef __INIT_INCLUDED__
#define __INIT_INCLUDED__

#include "data.h"

extern void read_parameters(const char *fname, params *pptr, connect *cptr);
extern void build_network(data *dptr, params *pptr, connect *cptr, cArray **combptr, cArray **hdptr);

extern void set_connectivity(connect *cptr, params *pptr);

extern void set_HD_COMB_weights(params *pptr, data *dptr, connect *cptr);

extern void set_COMB_HD_weights(params *pptr, data *dptr, connect *cptr);

extern void set_ROT_COMB_weights(params *pptr, data *dptr, connect *cptr);

extern void set_favoured_view(data *dptr, params *pptr);

extern void free_memory(params *pptr, data *dptr, connect *cptr, cArray **combptr, cArray **hdptr);

#endif