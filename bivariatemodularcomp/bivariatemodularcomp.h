#ifndef __BIVARIATEMODULARCOMP_H__
#define __BIVARIATEMODULARCOMP_H__
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_p.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include "mosaic_toeplitz_mul_ZZ_p.h"

class BivariateModularComp : public mosaic_toeplitz_mul_ZZ_p {
  NTL::ZZ_pX f_field;
  NTL::ZZ_pX F_field;
  NTL::Mat<NTL::ZZ_pX> B; // the right side of the matrix multiplication
  long sqrtP; // ceiling of the sqrt of the number of blocks
  

  /*******************************************************************
  * Helpers                                                          *
  *******************************************************************/
  
  // given the result of create_lhs_list, returns a matrix representation
  // with dimension ceil(len(type)/sqrtP) x sqrtP
  void create_lhs_matrix (NTL::Mat<NTL::ZZ_pX>&, const NTL::Vec<NTL::ZZ_pX>&);
  // creates a matrix of f^i for i = 0..sqrtP-1, and partitions each f^i
  // into len(type)-slices (number of blocks)
  void create_rhs_matrix (NTL::Mat<NTL::ZZ_pX>&, const NTL::Vec<NTL::ZZ_pX>&);
  // converts the sliced up functions into one
  void deslice (NTL::Vec<NTL::ZZ_pX>&, const NTL::Mat<NTL::ZZ_pX> &);

public:
  // the constructor takes in the polynomial (full precision), type and precision the answer should be
  BivariateModularComp(const NTL::ZZ_pX& f, const NTL::Vec<long> &type, long prec);


  // default ctor; does nothing
  BivariateModularComp();

  // initializes the function with type
  void init(const NTL::ZZ_pX& f, const NTL::Vec<long> &type, long prec);

   // multiplies rhs with the matrix created by the powers of f in
  // the given type by using Horner's rule
  NTL::Vec<NTL::ZZ_p> mult_right_Horners(const NTL::Vec<NTL::ZZ_p> &rhs);

  // multiply using the new algorithm
  NTL::Vec<NTL::ZZ_p> mult_right_comp (const NTL::Vec<NTL::ZZ_p> &rhs);
  
  // header for multiplying
  NTL::Vec<NTL::ZZ_p> mult_right (const NTL::Vec<NTL::ZZ_p> &rhs);
};
#endif
