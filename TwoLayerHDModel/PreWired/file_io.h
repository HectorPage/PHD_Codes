#ifndef __IO_INCLUDED__
#define __IO_INCLUDED__

#include "data.h"

extern FILE* file_open(const char *fname, const char *format);

extern void trim_space(char*str);

extern void save_rates_activations(data *dptr, params *pptr);

extern void save_pvector(data *dptr);

extern void save_initialiser(data *dptr, params *pptr);

#endif