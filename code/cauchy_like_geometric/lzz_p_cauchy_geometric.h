#ifndef __CAUCHY_GEOMETRIC_H
#define __CAUCHY_GEOMETRIC_H

#include <NTL/lzz_pX.h>
#include <NTL/mat_lzz_p.h>

#include "lzz_pX_CRT.h"
#include "lzz_p_toeplitz.h"

#define THRESHOLD 100

NTL_CLIENT

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/* Cauchy matrices on geometric progressions          */
/*----------------------------------------------------*/
/*----------------------------------------------------*/
class lzz_p_cauchy_geometric{
 private:
  long m, n;
 
 public:
  zz_p u1, v1, rho, sqrt_rho;
  Vec<zz_p> vec_toeplitz;
  lzz_p_toeplitz t;
  zz_pX_Multipoint_Geometric X, Y; // TODO: make it lzz_p_
  Vec<zz_p> powers_irho;

  /*---------------------------------------------------*/
  /* default constructor                               */
  /*---------------------------------------------------*/
  lzz_p_cauchy_geometric();

  /*---------------------------------------------------*/
  /* constructor                                       */
  /*---------------------------------------------------*/
  lzz_p_cauchy_geometric(const zz_p& a1, const zz_p& b1, const zz_p& rho, long mm, long nn);

  /*---------------------------------------------------*/
  /* dimensions                                        */
  /*---------------------------------------------------*/
  long NumRows() const;
  long NumCols() const;

  /*---------------------------------------------------*/
  /* computes output = M*input                         */
  /*---------------------------------------------------*/
  void mul_right(Vec<zz_p>& output, const Vec<zz_p>& input) const;
  void mul_right(Mat<zz_p>& output, const Mat<zz_p>& input) const;

  /*---------------------------------------------------*/
  /* computes output = M*input, without the diagonal   */
  /*---------------------------------------------------*/
  void mul_right_simple(Vec<zz_p>& output, const Vec<zz_p>& input) const;
  void mul_right_simple(Mat<zz_p>& output, const Mat<zz_p>& input) const;

  /*---------------------------------------------------*/
  /* computes output = M^t*input                       */
  /*---------------------------------------------------*/
  void mul_left(Vec<zz_p>& output, const Vec<zz_p>& input) const;
  void mul_left(Mat<zz_p>& output, const Mat<zz_p>& input) const;

  /*---------------------------------------------------*/
  /* computes output = M^t*input, without the diagonal */
  /*---------------------------------------------------*/
  void mul_left_simple(Vec<zz_p>& output, const Vec<zz_p>& input) const;
  void mul_left_simple(Mat<zz_p>& output, const Mat<zz_p>& input) const;

  /*---------------------------------------------------*/
  /* M as a dense matrix                               */
  /*---------------------------------------------------*/
  void to_dense(Mat<zz_p>& M) const;

  /*---------------------------------------------------*/
  /* X and Y are only built if needed                  */
  /*---------------------------------------------------*/
  void build_X_Y();

};


/*---------------------------------------------------*/
/* computes                                          */
/* 1/(u1-v1 rho^(-m+1)) ... 1/(u1-v1 rho^(n-1))      */
/* these are the entries of the toeplitz matrix      */
/* (with m rows and n columns)                       */
/*---------------------------------------------------*/
void prepare_inverses_cauchy(Vec<zz_p>& inverses, const zz_p& u1, const zz_p& v1, const zz_p& rho, long m, long n);



/*----------------------------------------------------*/
/*----------------------------------------------------*/
/* Cauchy-like matrices on geometric progressions     */
/*----------------------------------------------------*/
/*----------------------------------------------------*/

class lzz_p_cauchy_like_geometric{
public:
  lzz_p_cauchy_geometric C;
  Mat<zz_p> G, H;

  /*---------------------------------------------------*/
  /* default constructor                               */
  /*---------------------------------------------------*/
  lzz_p_cauchy_like_geometric();

  /*---------------------------------------------------*/
  /* constructor                                       */
  /*---------------------------------------------------*/
  lzz_p_cauchy_like_geometric(const Mat<zz_p>& U, const Mat<zz_p>& V, const zz_p& a1, const zz_p& b1, const zz_p& rho);

  /*---------------------------------------------------*/
  /* dimensions                                        */
  /*---------------------------------------------------*/
  long NumRows() const;
  long NumCols() const;
  long NumGens() const;

  /*---------------------------------------------------*/
  /* computes output = M*input                         */
  /*---------------------------------------------------*/
  void mul_right(Vec<zz_p>& output, const Vec<zz_p>& input) const;
  void mul_right_direct(Mat<zz_p> & out, const Mat<zz_p> & in) const ;
  void mul_right_sigma_UL(Mat<zz_p> & out, const Mat<zz_p> & in) ;
  void mul_right(Mat<zz_p>& output, const Mat<zz_p>& input) ;

  /*---------------------------------------------------*/
  /* computes output = M^t*input                       */
  /*---------------------------------------------------*/
  void mul_left(Vec<zz_p>& output, const Vec<zz_p>& input) const;
  void mul_left(Mat<zz_p>& output, const Mat<zz_p>& input) const;

  /*---------------------------------------------------*/
  /* M as a dense matrix                               */
  /*---------------------------------------------------*/
  void to_dense(Mat<zz_p>& M) const;
};

/*---------------------------------------------------*/
/* iCL = CL^(-1), represented for the dual operator  */
/* O(alpha n^2) algorithm                            */
/*---------------------------------------------------*/
long invert_direct(lzz_p_cauchy_like_geometric& Cinv, const lzz_p_cauchy_like_geometric& C);
long invert_block(lzz_p_cauchy_like_geometric& Cinv, const lzz_p_cauchy_like_geometric& C);
long invert(lzz_p_cauchy_like_geometric& Cinv, const lzz_p_cauchy_like_geometric& C);

long invert_fast(lzz_p_cauchy_like_geometric& Cinv, const lzz_p_cauchy_like_geometric& C, const long thresh = -1);


#endif
