#ifndef __hermite_pade_algebraic_h__
#define __hermite_pade_algebraic_h__

#include "bivariatemodularcomp.h"
#include <NTL/ZZ_pX.h>
#include <NTL/lzz_pX.h>
#include <NTL/tools.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include "mosaic_hankel.h"
#include "ZZ_pX_CRT.h"
#include "cauchy_geometric_special.h"
#include "hermite_pade.h"

class hermite_pade_algebraic : public hermite_pade{
  ZZX f_full_prec; // the polynomial in full precision
  ZZ denom = ZZ{1}; // the denominator for every entry in f

  /**** helpers *****************************************/ 
  
  void init(const ZZX &f, const Vec<long>& type, long prec, long fft_init);
 
  Vec<zz_pX> split_on_type(const Vec<zz_p> &v);
  
  void set_up_bmc() override;
  
  void init_bmc();
  
  public:
  /** Ctor ***********************************************
  * f: the generating series                             *
  * type: the type of the approximant                    *
  *******************************************************/
  hermite_pade_algebraic(const ZZX &f, const Vec<long>& type, long prec, long fft_init = 0);
  
  hermite_pade_algebraic(const ZZX &f, const ZZ &denom, const Vec<long>& type, long prec_inp, long fft_index);
  
  Vec<long> update_type();
  
  ~hermite_pade_algebraic(){}
};

#endif
