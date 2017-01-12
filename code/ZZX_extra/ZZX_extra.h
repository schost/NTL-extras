#ifndef ZZX_EXTRA__H
#define ZZX_EXTRA__H

#include <NTL/ZZX.h>

NTL_CLIENT

/*------------------------------------------------------------*/
/* random poly with deg(F,x) <= dx                            */
/* coefficients are in -2^{bound-1}..2^{bound-1}-1            */
/*------------------------------------------------------------*/
void random(ZZX & F, long bound, long dx);

/*------------------------------------------------------------*/
/* removes the contents of g and den                          */
/*------------------------------------------------------------*/
void simplify(ZZX & g, ZZ & den);

/*------------------------------------------------------------*/
/* finds ig, den_ig such that ig/den_ig * g/den_g = 1 mod x^t */
/*------------------------------------------------------------*/
void inverse_series(ZZX & ig, ZZ & den_ig, const ZZX & g, const ZZ & den_g, long t);

#endif
