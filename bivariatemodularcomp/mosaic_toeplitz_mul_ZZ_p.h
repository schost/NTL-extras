#ifndef __MOSAIC_TOEPLITZ_MUL_ZZ_P_H__
#define __MOSAIC_TOEPLITZ_MUL_ZZ_P_H__
#include <NTL/vector.h>
#include <NTL/ZZ_pX.h>

class mosaic_toeplitz_mul_ZZ_p{
protected:
  NTL::Vec<long> type; // the type for the matrix (in Hermite-Pade approximants)
  long prec; // precision of the returned vector
  bool initialized = false; 
  
  mosaic_toeplitz_mul_ZZ_p();
  mosaic_toeplitz_mul_ZZ_p(const NTL::Vec<long> &t,long p);
  
  // partitions the given vector into blocks from the given type and
  // converts each block into polynomials
  void create_lhs_list (NTL::Vec<NTL::ZZ_pX>&, const NTL::Vec<NTL::ZZ_p>&);

public:
	virtual NTL::Vec<NTL::ZZ_p> mult_right(const NTL::Vec<NTL::ZZ_p> &rhs) = 0;
	
	virtual ~mosaic_toeplitz_mul_ZZ_p(){}
};

#endif
