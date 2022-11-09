#ifndef __lagranzh__
#define __lagranzh__
//#include "fraction.h"
    
void mat_prod_r(Frac *a, Frac *b, int n);
void mat_prod_l(Frac *a, Frac *b, int n);
void transf(Frac *a, Frac *tr, int n);
void swap(Frac *a, int n, int p, int q, Frac *tr);

triple Lagranzh(Frac *a);

#endif