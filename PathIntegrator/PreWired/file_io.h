#ifndef __IO_INCLUDED__
#define __IO_INCLUDED__

#include "data.h"

extern FILE* file_open(const char *fname, const char *format);

extern void trim_space(char*str);

extern void save_rates_activations(data *dptr, params *pptr);

extern void save_input(data *dptr, params *pptr);

extern void save_RCweights_final(connect *cptr, params *pptr);

extern void save_RCweights_init(connect *cptr, params *pptr);

extern void load_full_RC(connect *cptr, params *pptr);

extern void save_RCweights_intermediate(connect *cptr, params *pptr);

extern void save_weight_vectors(data *dptr, params *pptr);

#endif