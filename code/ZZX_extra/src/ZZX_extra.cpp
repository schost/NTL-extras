#include <NTL/ZZX.h>
#include "ZZX_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* random poly with deg(F,x) <= dx                            */
/* coefficients are in -2^{bound-1}..2^{bound-1}-1            */
/*------------------------------------------------------------*/
void random(ZZX & F, long bound, long dx){
  if (dx == -1){
    clear(F);
    return;
  }

  ZZ half = ZZ{1} << (bound - 1);
  F.rep.SetLength(dx+1);
  for (long j = 0; j <= dx; j++)
    F.rep[j] = RandomBits_ZZ(bound) - half;
  F.normalize();
}

/*------------------------------------------------------------*/
/* removes the contents of g and den                          */
/*------------------------------------------------------------*/
void simplify(ZZX & g, ZZ & den){
  if (g == 0){
    den = 1;
    return;
  }
  ZZ h = GCD(den, coeff(g, 0));
  for (long i = 1; i <= deg(g); i++)
    h = GCD(h, coeff(g, i));
  den /= h;
  for (long i = 0; i <= deg(g); i++)
    SetCoeff(g, i, coeff(g, i) / h);
}

/*------------------------------------------------------------*/
/* finds ig, den_ig such that ig/den_ig * g/den_g = 1 mod x^t */
/*------------------------------------------------------------*/
void inverse_series(ZZX & ig, ZZ & den_ig, const ZZX & g, const ZZ & den_g, long t){
  clear(ig);
  SetCoeff(ig, 0, 1);
  den_ig = coeff(g, 0);
  long prec = 1;

  while (prec < t){
    prec = 2*prec;
    ZZX tmp = trunc(ig * g, prec);
    ZZX tmp2 = ZZX{2*den_ig} - tmp;
    ig = trunc(ig * tmp2, prec);
    den_ig = den_ig * den_ig;
    simplify(ig, den_ig);
  }

  ig = trunc(ig, t);
  ig *= den_g;
  simplify(ig, den_ig);
}
