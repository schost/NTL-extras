#ifndef __BIVARIATEMODULARCOMP_H__
#define __BIVARIATEMODULARCOMP_H__
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_p.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>

class BivariateModularComp {
  NTL::ZZ_pX f_field;
  NTL::ZZ_pX F_field;
  NTL::Vec<NTL::ZZ_pX> fs_field;
  NTL::Vec<long> type; // the type for the matrix (in Hermite-Pade approximants)
  NTL::Mat<NTL::ZZ_pX> B; // the right side of the matrix multiplication
  long prec; // precision of the returned vector
  long sqrtP; // ceiling of the sqrt of the number of blocks
  bool initialized = false;
  /* switches which multiplication to use: */
  /* 0 - bivariate modular composition     */
  /* 1 - naive multiplication              */
  long mode; 

  /*******************************************************************
  * Helpers                                                          *
  *******************************************************************/
  // partitions the given vector into blocks from the given type and
  // converts each block into polynomials
  void create_lhs_list (NTL::Vec<NTL::ZZ_pX>&, const NTL::Vec<NTL::ZZ_p>&);
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

  // construct given sufficient powers of f, see init
  BivariateModularComp(const NTL::Vec<NTL::ZZ_pX> &fs, const NTL::Vec<long> &type, long prec);

  // default ctor; does nothing
  BivariateModularComp();

  // initializes the function with type
  void init(const NTL::ZZ_pX& f, const NTL::Vec<long> &type, long prec);

  // initialize without calculating the powers, assumes fs = {f^0, f^1, f^2 ...}
  // throws an exception if not enough powers of f are given
  void init (const NTL::Vec<NTL::ZZ_pX>& fs, const NTL::Vec<long> &type, long prec);
  
   // multiplies rhs with the matrix created by the powers of f in
  // the given type by using Horner's rule
  NTL::Vec<NTL::ZZ_p> mult_right_Horners(const NTL::Vec<NTL::ZZ_p> &rhs);

  // multiply using the new algorithm
  NTL::Vec<NTL::ZZ_p> mult_right_comp (const NTL::Vec<NTL::ZZ_p> &rhs);
  
  // multiply naively as a0*f0 + a1*f1 + a2*f2
  NTL::Vec<NTL::ZZ_p> mult_right_naive (const NTL::Vec<NTL::ZZ_p> &rhs);
  
  // header for multiplying
  NTL::Vec<NTL::ZZ_p> mult_right (const NTL::Vec<NTL::ZZ_p> &rhs);
};
#endif
