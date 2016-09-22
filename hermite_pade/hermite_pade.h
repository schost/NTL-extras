#ifndef __hermite_pade_h__
#define __hermite_pade_h__

#include "bivariatemodularcomp.h"
#include <NTL/ZZ_pX.h>
#include <NTL/zz_pX.h>
#include <NTL/tool.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>

class hermite_pade{
  bivariatemodularcomp M; // the actually hermite-pade matrix
  Vec<zz_p> e, f; // for precondition the Cauchy Matrix
  Vec<zz_pX_Multipoint_Geometric> vec_X_int, vec_Y_int; // to store already calculated powers

  public:
  /** Ctor ***********************************************
  * f: the generating series                             *
  * type: the type of the approximant                    *
  *******************************************************/
  hermite_pade(const ZZX &f, const Vec<long>& type){
    M.init(f,type,deg(f)+1);
    // create the mosaic hankel matrices
    // convert to cauchy matrix
  }

  void get_rand_sol(Vec<ZZ> &);
}

#endif
