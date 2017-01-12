#ifndef __LZZ_P_TOEPLITZ_LIKE_H
#define __LZZ_P_TOEPLITZ_LIKE_H

#include <NTL/lzz_pX.h>
#include <NTL/mat_lzz_p.h>

NTL_CLIENT

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/* toeplitz-like matrices, for the operator phi_minus */
/* (cf. Kaltofen 1994)                                */
/*----------------------------------------------------*/
/*----------------------------------------------------*/

class lzz_p_toeplitz_like{

 public:
  Mat<zz_p> G, H;
  long alpha, m, n; // generators of size (m x alpha) and (n x alpha)
  
  lzz_p_toeplitz_like(const Mat<zz_p>& G0, const Mat<zz_p>& H0);

  /*----------------------------------------------------*/
  /* output = M^t * input                               */
  /*----------------------------------------------------*/
  void mul_left(Vec<zz_p>& output, const Vec<zz_p>& input);

  /*----------------------------------------------------*/
  /* output = M * input                                 */
  /*----------------------------------------------------*/
  void mul_right(Vec<zz_p>& output, const Vec<zz_p>& input);

  /*----------------------------------------------------*/
  /* reconstructs M from its generators                 */
  /*----------------------------------------------------*/
  void to_dense(Mat<zz_p>& output);


  void mul_right_dac(Mat<zz_p>& output, const Mat<zz_p>& input);
};

/*----------------------------------------------------*/
/* down-shift matrix of size n                        */
/*----------------------------------------------------*/
Mat<zz_p> do_Z(long n);

/*----------------------------------------------------*/
/* M - Z^t M Z = M - (M shifted up, left)             */
/*----------------------------------------------------*/
void phi_minus(Mat<zz_p>& output, const Mat<zz_p>& input);





#endif
