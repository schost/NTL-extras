#ifndef __NON_ALGEBRAIC_TOEPLITZ_MUL_ZZ_P_H_
#define __NON_ALGEBRAIC_TOEPLITZ_MUL_ZZ_P_H_
#include "mosaic_toeplitz_mul_ZZ_p.h"

class non_algebraic_mosaic_toeplitz_mul_ZZ_p: public mosaic_toeplitz_mul_ZZ_p{
	NTL::Vec<NTL::ZZ_pX> fs_field;
public:
	void init (const NTL::Vec<NTL::ZZ_pX>& fs, const NTL::Vec<long> &type, long prec);

  non_algebraic_mosaic_toeplitz_mul_ZZ_p(const NTL::Vec<NTL::ZZ_pX> &fs, const NTL::Vec<long> &type, long prec);
  
  NTL::Vec<NTL::ZZ_p> mult_right(const NTL::Vec<NTL::ZZ_p> &rhs);
};

#endif
