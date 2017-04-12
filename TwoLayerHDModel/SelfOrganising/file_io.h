#ifndef __IO_INCLUDED__
#define __IO_INCLUDED__

#include "data.h"

extern FILE* file_open(const char *fname, const char *format);

extern void trim_space(char*str);

extern void save_rates_activations(data *dptr, params *pptr);

extern void save_initialiser(data *dptr, params *pptr);

extern void save_HDtoCOMB_weights(connect *cptr, params *pptr);

extern void save_COMBtoHD_weights(connect *cptr, params *pptr);

extern void load_full_HDtoCOMB(connect *cptr, params *pptr);

extern void save_ROTtoCOMB_weights(connect *cptr, params *pptr);

extern void save_initial_HDtoCOMB_weights(connect *cptr, params *pptr);

extern void save_initial_COMBtoHD_weights(connect *cptr, params *pptr);

extern void save_initial_ROTtoCOMB_weights(connect *cptr, params *pptr);


#endif