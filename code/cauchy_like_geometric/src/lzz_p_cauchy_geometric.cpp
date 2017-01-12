#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>

#include "vec_lzz_p_extra.h"
#include "lzz_pX_CRT.h"
#include "lzz_p_toeplitz.h"
#include "lzz_p_cauchy_geometric.h"

NTL_CLIENT

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/* Cauchy matrices on geometric progressions          */
/*----------------------------------------------------*/
/*----------------------------------------------------*/

/*---------------------------------------------------*/
/* default constructor                               */
/*---------------------------------------------------*/
lzz_p_cauchy_geometric::lzz_p_cauchy_geometric() {
  m = 0;
  n = 0;
}

/*---------------------------------------------------*/
/* constructor                                       */
/*---------------------------------------------------*/
lzz_p_cauchy_geometric::lzz_p_cauchy_geometric(const zz_p& a1, const zz_p& b1, const zz_p& q, long mm, long nn){
  u1 = a1;
  v1 = b1;
  rho = q;
  sqrt_rho = 0;
  m = mm;
  n = nn;

  if (m == 0 && n == 0){
    t = lzz_p_toeplitz();
    powers_irho.SetLength(0);
  }
  else{
    prepare_inverses_cauchy(vec_toeplitz, u1, v1, rho, m, n); 
   inverse_powers(powers_irho, rho, m);
    t = lzz_p_toeplitz(vec_toeplitz, m, n);
  }
}

/*---------------------------------------------------*/
/* dimensions                                        */
/*---------------------------------------------------*/
long lzz_p_cauchy_geometric::NumRows() const {
  return m;
}

long lzz_p_cauchy_geometric::NumCols() const {
  return n;
}

/*---------------------------------------------------*/
/* computes output = M*input                         */
/*---------------------------------------------------*/
void lzz_p_cauchy_geometric::mul_right(Vec<zz_p>& output, const Vec<zz_p>& input) const {
  t.mul_right(output, input);
  long m = NumRows();
  for (long i = 0; i < m; i++)
    output[i] *= powers_irho[i];
}

/*---------------------------------------------------*/
/* computes output = M*input                         */
/*---------------------------------------------------*/
void lzz_p_cauchy_geometric::mul_right(Mat<zz_p>& output, const Mat<zz_p>& input) const {
  t.mul_right(output, input);
  long m = NumRows();
  long a = input.NumCols();
  for (long i = 0; i < m; i++)
    for (long j = 0; j < a; j++)
      output[i][j] *= powers_irho[i];
}

/*---------------------------------------------------*/
/* computes output = M*input without the diagonal    */
/*---------------------------------------------------*/
void lzz_p_cauchy_geometric::mul_right_simple(Vec<zz_p>& output, const Vec<zz_p>& input) const {
  t.mul_right(output, input);
}

/*---------------------------------------------------*/
/* computes output = M*input without the diagonal    */
/*---------------------------------------------------*/
void lzz_p_cauchy_geometric::mul_right_simple(Mat<zz_p>& output, const Mat<zz_p>& input) const {
  t.mul_right(output, input);
}

/*---------------------------------------------------*/
/* computes output = M^t*input                       */
/*---------------------------------------------------*/
void lzz_p_cauchy_geometric::mul_left(Vec<zz_p>& output, const Vec<zz_p>& input) const {
  Vec<zz_p> new_in = input;
  long m = NumRows();
  for (long i = 0; i < m; i++)
    new_in[i] *= powers_irho[i];
  t.mul_left(output, new_in);
}

/*---------------------------------------------------*/
/* computes output = M^t*input                       */
/*---------------------------------------------------*/
void lzz_p_cauchy_geometric::mul_left(Mat<zz_p>& output, const Mat<zz_p>& input) const {
  Mat<zz_p> new_in = input;
  long m = NumRows();
  long a = input.NumCols();
  for (long i = 0; i < m; i++){
    zz_p *elts = new_in[i].elts();
    for (long j = 0; j < a; j++)
      elts[j] *= powers_irho[i];
  }
  t.mul_left(output, new_in);
}

/*---------------------------------------------------*/
/* computes output = M^t*input without the diagonal  */
/*---------------------------------------------------*/
void lzz_p_cauchy_geometric::mul_left_simple(Vec<zz_p>& output, const Vec<zz_p>& input) const {
  t.mul_left(output, input);
}

/*---------------------------------------------------*/
/* computes output = M^t*input without the diagonal  */
/*---------------------------------------------------*/
void lzz_p_cauchy_geometric::mul_left_simple(Mat<zz_p>& output, const Mat<zz_p>& input) const {
  t.mul_left(output, input);
}

/*---------------------------------------------------*/
/* M as a dense matrix                               */
/*---------------------------------------------------*/
void lzz_p_cauchy_geometric::to_dense(Mat<zz_p>& M) const {
  long m = NumRows();
  long n = NumCols();
  M.SetDims(m, n);
  for (long i = 0; i < m; i++)
    for (long j = 0; j < n; j++)
      M[i][j] = to_zz_p(1) / (u1*power(rho, i)-v1*power(rho, j));
}

/*----------------------------------------------------*/
/* builds the left / right vandermonde matrices       */
/* we do it only if they are needed                   */
/*----------------------------------------------------*/
void lzz_p_cauchy_geometric::build_X_Y(){
  if (X.n != 0)
    return;

  // TODO: check if it exists!
  if (sqrt_rho == 0){
    long rho_l = rho._zz_p__rep;
    sqrt_rho = to_zz_p(to_long(SqrRootMod(to_ZZ(rho_l), to_ZZ(zz_p::modulus()))));
  }

  // slightly faster in this case
  if (u1 == to_zz_p(1))
    X = zz_pX_Multipoint_Geometric(sqrt_rho, m, 0);
  else
    X = zz_pX_Multipoint_Geometric(sqrt_rho, m, u1);

  // slightly faster in this case
  if (v1 == to_zz_p(1))
    Y = zz_pX_Multipoint_Geometric(sqrt_rho, n, 0);
  else
    Y = zz_pX_Multipoint_Geometric(sqrt_rho, n, v1);
}

/*---------------------------------------------------*/
/* computes                                          */
/* 1/(u1-v1 rho^(-m+1)) ... 1/(u1-v1 rho^(n-1))      */
/* these are the entries of the toeplitz matrix      */
/* (with m rows and n columns)                       */
/*---------------------------------------------------*/
void prepare_inverses_cauchy(Vec<zz_p>& inverses, const zz_p& u1, const zz_p& v1, const zz_p& rho, long m, long n){
  
  Vec<zz_p> vec_den;
  zz_p irho = 1/rho;
  vec_den.SetLength(m+n-1);
  if (m+n-1 == 0)
    return;
  vec_den[0] = -v1*power(irho, m-1);
  for (long i = 1; i < m+n-1; i++){
    vec_den[i] = vec_den[i-1]*rho;
    vec_den[i-1] += u1;
  }
  vec_den[n+m-2] += u1;
  inv(inverses, vec_den);
}
