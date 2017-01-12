#ifndef __ZZ_P_BLOCK_SYLVESTER_H__
#define __ZZ_P_BLOCK_SYLVESTER_H__

#include <NTL/ZZX.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>

NTL_CLIENT

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/* matrices of the form [S0 .. St]                    */
/* each Si is Sylvester for some polynomial fi        */
/*----------------------------------------------------*/
/*----------------------------------------------------*/
class ZZ_p_block_sylvester{ 
protected:
  Vec<long> type; // the type for the matrix (in Hermite-Pade approximants)
  long max_of_type;
  long prec; // precision of the returned vector
  bool initialized = false; 
  
  ZZ_p_block_sylvester();
  ZZ_p_block_sylvester(const Vec<long> &t, long prec);

  /*----------------------------------------------------*/
  /* partitions the given vector into blocks            */
  /* converts each block into polynomials               */
  /*----------------------------------------------------*/
  void create_lhs_list (Vec<ZZ_pX>&, const Vec<ZZ_p>&);

public:
  virtual Vec<ZZ_p> mul_right(const Vec<ZZ_p> &rhs) = 0;
  virtual ~ZZ_p_block_sylvester(){}
};


/*----------------------------------------------------*/
/*----------------------------------------------------*/
/* the case of general fi's                           */
/*----------------------------------------------------*/
/*----------------------------------------------------*/
class ZZ_p_block_sylvester_general: public ZZ_p_block_sylvester {
  Vec<ZZ_pX> f;

public:
  /*----------------------------------------------------*/
  /* is this useful?                                    */
  /*----------------------------------------------------*/
  void init (const Vec<ZZ_pX> &fs, const Vec<long> &type, long prec);

  /*----------------------------------------------------*/
  /* input: Vec of polynomials fs                       */
  /*        type                                        */
  /*        output precision                            */
  /*----------------------------------------------------*/
  ZZ_p_block_sylvester_general(const Vec<ZZ_pX> &fs, const Vec<long> &type, long prec);

  /*----------------------------------------------------*/
  /* right multiplication                               */
  /*----------------------------------------------------*/
  Vec<ZZ_p> mul_right(const Vec<ZZ_p> &rhs); // TODO: mult -> mul  // TODO: void version
};






/*----------------------------------------------------*/
/*----------------------------------------------------*/
/* the case with fi = f^i                             */
/*----------------------------------------------------*/
/*----------------------------------------------------*/
class ZZ_p_bivariate_modular_composition : public ZZ_p_block_sylvester {

  ZZ_pX f_field;
  ZZ_pX F_field;
  Mat<ZZ_pX> B; // the right side of the matrix multiplication
  long sqrtP; // ceiling of the sqrt of the number of blocks
  
  /*----------------------------------------------------*/
  /* given the result of create_lhs_list, returns a     */
  /* matrix representation with dimension               */
  /*  ceil(len(type)/sqrtP) x sqrtP                     */
  /*----------------------------------------------------*/
  void create_lhs_matrix (Mat<ZZ_pX>& mat, const Vec<ZZ_pX>& vec);

  /*----------------------------------------------------*/
  /* creates a matrix of f^i for i = 0..sqrtP-1         */
  /* partitions each f^i into len(type) slices          */
  /*----------------------------------------------------*/
  void create_rhs_matrix (Mat<ZZ_pX>&, const Vec<ZZ_pX>&);

  /*----------------------------------------------------*/
  /* converts the sliced up matrix into a vector        */
  /*----------------------------------------------------*/
  void deslice (Vec<ZZ_pX>& de_sliced, const Mat<ZZ_pX> & sliced_mat);

public:
  /*----------------------------------------------------*/
  /* ?                                                  */
  /*----------------------------------------------------*/
  void init(const ZZ_pX& f, const Vec<long> &type, long prec);

  /*----------------------------------------------------*/
  /* input: polynomial f                                */
  /*        type                                        */
  /*        output precision                            */
  /*----------------------------------------------------*/
  ZZ_p_bivariate_modular_composition(const ZZ_pX& f, const Vec<long> &type, long prec);

  /*----------------------------------------------------*/
  /* default constructor; does nothing                  */
  /*----------------------------------------------------*/
  ZZ_p_bivariate_modular_composition();

  /*----------------------------------------------------*/
  /* multiplies rhs using Horner's rule                 */
  /*----------------------------------------------------*/
  Vec<ZZ_p> mul_right_Horners(const Vec<ZZ_p> &rhs);  //TODO: void version

  /*----------------------------------------------------*/
  /* using the baby steps / giant steps algorithm       */
  /*----------------------------------------------------*/
  Vec<ZZ_p> mul_right_comp (const Vec<ZZ_p> &rhs);  //TODO: void version
  
  /*----------------------------------------------------*/
  /* header for multiplying                             */
  /*----------------------------------------------------*/
  Vec<ZZ_p> mul_right (const Vec<ZZ_p> &rhs);  //TODO: void version
};


#endif
