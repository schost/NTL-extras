#ifndef __BIVARIATEMODULARCOMP_H__
#define __BIVARIATEMODULARCOMP_H__
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_p.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>

class BivariateModularComp {
  NTL::ZZX f; // the polynomial
  NTL::Vec<NTL::ZZX> fs; // the list for f^i for i = 0..sqrtP-1
  NTL::ZZX F; // the largest power of f
  NTL::ZZ_pX f_field;
  NTL::ZZ_pX F_field;
  NTL::Vec<NTL::ZZ_pX> fs_field;
  NTL::Vec<long> type; // the type for the matrix (in Hermite-Pade approximants)
  NTL::Mat<NTL::ZZ_pX> B; // the right side of the matrix multiplication
  unsigned long prec; // precision of the returned vector
  unsigned long sqrtP; // ceiling of the sqrt of the number of blocks
  unsigned long totalVars = 0; // total number of variables
  bool initialized = false;

  /*******************************************************************
  * Helpers                                                          *
  *******************************************************************/
  // partitions the given vector into blocks from the given type and
  // converts each block into polynomials
  NTL::Vec<NTL::ZZ_pX> create_lhs_list (const NTL::Vec<NTL::ZZ_p>&);
  // given the result of create_lhs_list, returns a matrix representation
  // with dimension ceil(len(type)/sqrtP) x sqrtP
  void create_lhs_matrix (NTL::Mat<NTL::ZZ_pX>&,const NTL::Vec<NTL::ZZ_pX>&);
  // creates a matrix of f^i for i = 0..sqrtP-1, and partitions each f^i
  // into len(type)-slices (number of blocks)
  void create_rhs_matrix (NTL::Mat<NTL::ZZ_pX>&,const NTL::Vec<NTL::ZZ_pX>&);
  // converts the sliced up functions into one
  void deslice (NTL::Mat<NTL::ZZ_pX>&,const NTL::Mat<NTL::ZZ_pX> &);

public:
  // the constructor takes in the polynomial (full precision), type and precision the answer should be
  BivariateModularComp(const NTL::ZZX& f, const NTL::Vec<long> &type, unsigned long prec);

  // default ctor; does nothing
  BivariateModularComp();

  // initializes the function with type
  void init(const NTL::ZZX& f, const NTL::Vec<long> &type, unsigned long prec);

  // multiplies rhs with the matrix created by the powers of f in
  // the given type by using Horner's rule
  NTL::Vec<NTL::ZZ_p> mult_Horners(const NTL::Vec<NTL::ZZ_p> &rhs);

  // multiply using the new algorithm
  NTL::Vec<NTL::ZZ_p> mult (const NTL::Vec<NTL::ZZ_p> &rhs);
};

// return f^p
NTL::ZZ_pX pow (const NTL::ZZ_pX& f, unsigned long p);
#endif
