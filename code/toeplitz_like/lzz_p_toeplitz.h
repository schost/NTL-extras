#ifndef __LZZ_P_TOEPLITZ_H
#define __LZZ_P_TOEPLITZ_H

#include <NTL/lzz_pX.h>
#include <NTL/mat_lzz_p.h>

#include "lzz_pX_CRT.h"


NTL_CLIENT

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/* Toeplitz matrices                                  */
/* stored as                                          */
/*       a2 a3 a4 a5                                  */
/*       a1 a2 a3 a4                                  */
/*       a0 a1 a2 a3                                  */
/*----------------------------------------------------*/
/*----------------------------------------------------*/

class lzz_p_toeplitz{
 public:
  // n rows, m columns
  long n, m;
  Vec<zz_p> data;
  Vec<zz_p> data_rev;
  fftRep fft_data;

  zz_pX_Multipoint_CTFT c;
  Vec<zz_p> TFT_data;

  /*----------------------------------------------------*/
  /* default constructor                                */
  /*----------------------------------------------------*/
  lzz_p_toeplitz();

  /*----------------------------------------------------*/
  /* input vector is as showed above                    */
  /*----------------------------------------------------*/
  lzz_p_toeplitz(const Vec<zz_p>& input, long rows, long cols);

  /*---------------------------------------------------*/
  /* dimensions                                        */
  /*---------------------------------------------------*/
  long NumRows() const;
  long NumCols() const;

  /*---------------------------------------------------*/
  /* data access                                       */
  /*---------------------------------------------------*/
  const zz_p& operator ()(long i, long j) const;

  /*---------------------------------------------------*/
  /* computes output = M*input                         */
  /*---------------------------------------------------*/
  void mul_right(Vec<zz_p>& output, const Vec<zz_p>& input) const;
  void mul_right(Mat<zz_p>& output, const Mat<zz_p>& input) const;

  void mul_right_CTFT(Vec<zz_p>& res, const Vec<zz_p>& input) const;

  /*---------------------------------------------------*/
  /* computes output = M^t*input                       */
  /*---------------------------------------------------*/
  void mul_left(Vec<zz_p>& output, const Vec<zz_p>& input) const;
  void mul_left(Mat<zz_p>& output, const Mat<zz_p>& input) const;

  /*---------------------------------------------------*/
  /* M as a dense matrix                               */
  /*---------------------------------------------------*/
  void to_dense(Mat<zz_p>& Mdense) const;
};

#endif
