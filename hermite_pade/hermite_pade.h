#ifndef __hermite_pade_h__
#define __hermite_pade_h__

#include "bivariatemodularcomp.h"
#include <NTL/ZZ_pX.h>
#include <NTL/lzz_pX.h>
#include <NTL/tools.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include "mosaic_hankel.h"
#include "ZZ_pX_CRT.h"
#include "cauchy_geometric_special.h"

class hermite_pade{
  Vec<BivariateModularComp> vec_M;
  cauchy_like_geometric_special invM; //the inverse mod p
  Vec<ZZ_p> e, f; // for precondition the Cauchy Matrix
  Vec<ZZ_pX_Multipoint_FFT> vec_X_int, vec_Y_int; // to store already calculated powers
  BivariateModularComp *M; // the current M
  ZZ_pX_Multipoint_FFT *X_int, *Y_int; // the current preconditioner
  ZZX f_full_prec; // the polynomial in full precision
  Vec<long> type; // the type of the approximant
  long rank; // the rank of M

  /**** helpers *****************************************/ 

  public:
  /** Ctor ***********************************************
  * f: the generating series                             *
  * type: the type of the approximant                    *
  *******************************************************/
  hermite_pade(const ZZX &f, const Vec<long>& type);

  // computes a random solution to the system
  void get_rand_sol(Vec<ZZ> &);
};

#endif
