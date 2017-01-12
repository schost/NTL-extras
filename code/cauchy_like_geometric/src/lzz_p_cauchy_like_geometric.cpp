#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>
#include <assert.h>

#include "vec_lzz_p_extra.h"
#include "lzz_pX_CRT.h"
#include "lzz_p_toeplitz.h"
#include "lzz_p_toeplitz_like.h"
#include "lzz_p_cauchy_geometric.h"

NTL_CLIENT

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/* Cauchy-like matrices on geometric progressions     */
/*----------------------------------------------------*/
/*----------------------------------------------------*/

/*---------------------------------------------------*/
/* default constructor                               */
/*---------------------------------------------------*/
lzz_p_cauchy_like_geometric::lzz_p_cauchy_like_geometric(){ 
}
  
/*---------------------------------------------------*/
/* constructor                                       */
/*---------------------------------------------------*/
lzz_p_cauchy_like_geometric::lzz_p_cauchy_like_geometric(const Mat<zz_p>& U, const Mat<zz_p>& V, 
							 const zz_p& a1, const zz_p& b1, const zz_p& rho){
  G = U;
  H = V;
  C = lzz_p_cauchy_geometric(a1, b1, rho, G.NumRows(), H.NumRows());
}

/*---------------------------------------------------*/
/* dimensions                                        */
/*---------------------------------------------------*/
long lzz_p_cauchy_like_geometric::NumRows() const {
  return G.NumRows();
}

long lzz_p_cauchy_like_geometric::NumCols() const {
  return H.NumRows();
}

long lzz_p_cauchy_like_geometric::NumGens() const {
  return G.NumCols();
}

/*---------------------------------------------------*/
/* computes output = M*input                         */
/*---------------------------------------------------*/
void lzz_p_cauchy_like_geometric::mul_right(Vec<zz_p>& output, const Vec<zz_p>& input) const {
  long m = NumRows();
  long n = NumCols();
  long alpha = NumGens();

  Vec<zz_p> new_in, new_out;
  new_in.SetLength(n);

  output.SetLength(m);
  for (long i = 0; i < m; i++)
    output[i] = 0;

  for (long i = 0; i < alpha; i++){
    for (long j = 0; j < n; j++)
      new_in[j] = H[j][i]*input[j];
    C.mul_right(new_out, new_in);
    for (long j = 0; j < m; j++)
      output[j] += G[j][i]*new_out[j];
  }
}

/*---------------------------------------------------*/
/* computes output = M*input                         */
/*---------------------------------------------------*/
void lzz_p_cauchy_like_geometric::mul_right_direct(Mat<zz_p>& output, const Mat<zz_p>& input) const {
  long m = NumRows();
  long n = NumCols();
  long alpha = NumGens();
  long beta = input.NumCols();
  output.SetDims(m, beta);

  if (m == 0 || n == 0){
    for (long i = 0; i < m; i++){
      zz_p * elts = output[i].elts();
      for (long j = 0; j < beta; j++)
	elts[j] = 0;
    }
    return;
  }

  Vec<zz_p> new_in, new_out, vec_out;
  new_in.SetLength(n);
  vec_out.SetLength(m);

  const lzz_p_toeplitz * t = &(C.t);

  for (long i = 0; i < beta; i++){
    for (long j = 0; j < m; j++)
      vec_out[j] = 0;
    for (long k = 0; k < alpha; k++){
      for (long j = 0; j < n; j++)
	new_in[j] = H[j][k]*input[j][i];
      t->mul_right(new_out, new_in);
      for (long j = 0; j < m; j++)
	vec_out[j] += G[j][k]*new_out[j];
    }
    const zz_p * diag = C.powers_irho.elts();
    for (long j = 0; j < m; j++)
      output[j][i] = diag[j]*vec_out[j];
  }
}

/*---------------------------------------------------*/
/* computes output = M*input                         */
/*---------------------------------------------------*/
void lzz_p_cauchy_like_geometric::mul_right_sigma_UL(Mat<zz_p> & output, const Mat<zz_p> & input) {

  long m = NumRows();
  long n = NumCols();
  long alpha = NumGens();
  long beta = input.NumCols();

  if (m == 0 || n == 0){
    for (long i = 0; i < m; i++){
      zz_p * elts = output[i].elts();
      for (long j = 0; j < beta; j++)
	elts[j] = 0;
    }
    return;
  }

  C.build_X_Y();

  Vec<zz_p> u;
  u.SetLength(m);
  u[0] = 1/C.u1;
  zz_p inv_rho = 1/C.rho;
  for (long i = 1; i < m; i++)
    u[i] = inv_rho*u[i-1];

  Vec<zz_p> v;
  v.SetLength(n);
  v[0] = power(C.v1, n);
  zz_p rho_n = power(C.rho, n);
  for (long i = 1; i < n; i++)
    v[i] = rho_n*v[i-1];

  Mat<zz_p> Gt, Ht;
  Gt.SetDims(m, alpha+2);
  Ht.SetDims(n, alpha+2);
  Vec<zz_p> in, out;
  in.SetLength(m);

  for (long i = 0; i < alpha; i++){
    for (long j = 0; j < m; j++)
      in[j] = u[j]*G[j][i];
    C.X.mul_left(out, in);
    for (long j = 0; j < m; j++)
      Gt[m-1-j][i] = out[j];
  }

  Vec<zz_p> out2;
  mul_right(out, v);
  for (long j = 0; j < m; j++)
    out[j] *= u[j];
  C.X.mul_left(out2, out);
  for (long j = 0; j < m; j++)
    Gt[m-1-j][alpha] = out2[j];

  for (long j = 0; j < m-1; j++)
    Gt[j][alpha+1] = 0;
  Gt[m-1][alpha+1] = 1;

  in.SetLength(n);
  for (long i = 0; i < alpha; i++){
    for (long j = 0; j < n; j++)
      in[j] = H[j][i];
    C.Y.mul_left(out, in);
    for (long j = 0; j < n; j++)
      Ht[j][i] = out[j];
  }

  for (long j = 0; j < n-1; j++)
    Ht[j][alpha] = 0;
  Ht[n-1][alpha] = 1;

  mul_left(out, u);
  C.Y.mul_left(out2, out);
  for (long j = 0; j < n-1; j++)
    Ht[j][alpha+1] = out2[j+1];
  Ht[n-1][alpha+1] = 0;

  lzz_p_toeplitz_like M(Gt, Ht);

  Mat<zz_p> M_input, M_output;
  M_input.SetDims(n, beta);

  in.SetLength(n);
  for (long i = 0; i < beta; i++){
    for (long j = 0; j < n; j++)
      in[j] = input[j][i];
    C.Y.inverse_mul_right(out, in);
    for (long j = 0; j < n; j++)
      M_input[j][i] = out[j];
  }

  M.mul_right_dac(M_output, M_input);

  output.SetDims(m, beta);
  in.SetLength(m);
  for (long i = 0; i < beta; i++){
    for (long j = 0; j < m; j++)
      in[j] = M_output[m-1-j][i];
    C.X.inverse_mul_left(out, in);
    for (long j = 0; j < m; j++)
      output[j][i] = out[j];
  }
}

/*---------------------------------------------------*/
/* computes output = M*input                         */
/*---------------------------------------------------*/
void lzz_p_cauchy_like_geometric::mul_right(Mat<zz_p> & output, const Mat<zz_p> & input) {
  long alpha = NumGens();
  if (alpha < 36)
    mul_right_direct(output, input);
  else
    mul_right_sigma_UL(output, input);
}


/*---------------------------------------------------*/
/* computes output = M^t*input                       */
/*---------------------------------------------------*/
void lzz_p_cauchy_like_geometric::mul_left(Vec<zz_p>& output, const Vec<zz_p>& input) const {

  long m = NumRows();
  long n = NumCols();
  long alpha = NumGens();

  Vec<zz_p> new_in, new_out;
  new_in.SetLength(m);

  output.SetLength(n);
  for (long i = 0; i < n; i++)
    output[i] = 0;

  for (long i = 0; i < alpha; i++){
    for (long j = 0; j < m; j++)
      new_in[j] = G[j][i]*input[j];
    C.mul_left(new_out, new_in);
    for (long j = 0; j < n; j++)
      output[j] += H[j][i]*new_out[j];
  }
}


/*---------------------------------------------------*/
/* computes output = M^t*input                       */
/*---------------------------------------------------*/
void lzz_p_cauchy_like_geometric::mul_left(Mat<zz_p>& output, const Mat<zz_p>& input) const {


  long m = NumRows();
  long n = NumCols();
  long alpha = NumGens();
  long beta = input.NumCols();
  output.SetDims(n, beta);

  if (m == 0 || n == 0){
    for (long i = 0; i < n; i++){
      zz_p * elts = output[i].elts();
      for (long j = 0; j < beta; j++)
	elts[j] = 0;
    }
    return;
  }

  Vec<zz_p> new_in, new_in_2, new_out, vec_out;
  new_in.SetLength(m);
  new_in_2.SetLength(m);
  vec_out.SetLength(n);

  for (long i = 0; i < beta; i++){
    const zz_p * diag = C.powers_irho.elts();
    for (long j = 0; j < m; j++)
      new_in_2[j] = input[j][i] * diag[j];
    for (long j = 0; j < n; j++)
      vec_out[j] = 0;
    for (long k = 0; k < alpha; k++){
      for (long j = 0; j < m; j++)
	new_in[j] = G[j][k]*new_in_2[j];
      C.mul_left_simple(new_out, new_in);
      for (long j = 0; j < n; j++)
	vec_out[j] += H[j][k]*new_out[j];
    }
    for (long j = 0; j < n; j++)
      output[j][i] = vec_out[j];
  }
}

/*---------------------------------------------------*/
/* M as a dense matrix                               */
/*---------------------------------------------------*/
void lzz_p_cauchy_like_geometric::to_dense(Mat<zz_p>& M) const {
  long m = NumRows();
  long n = NumCols();

  M = G*transpose(H);
  for (long i = 0; i < m; i++)
    for (long j = 0; j < n; j++)
      M[i][j] /= (C.u1*power(C.rho, i)-C.v1*power(C.rho, j));
}

/*---------------------------------------------------*/
/* needs much more data!                             */
/*---------------------------------------------------*/
static 
long threshold(long alpha){
  if (zz_p::modulus() < (1L << 20)){
    if (alpha < 10)
      return 200;
    if (alpha < 15)
      return 400;
    return 650;
  }
  else{
    if (alpha < 10)
      return 200;
    if (alpha < 20)
      return 400;
    return 700;

  }
}


/*------------------------------------------------------------*/
/* inverts a special-Cauchy-like matrix                       */
/* assumes generic rank profile                               */
/* the generators are Yp, Zp                                  */
/* alpha = displacement rank                                  */
/* overwrites Yp, Zp                                          */
/* i0, j0 are the initial indices in i, j resp.               */
/* i1, j1 are the end indices in i, j resp.                   */  
/*------------------------------------------------------------*/
static long invert_block_raw(Mat<zz_p>& Yp_out, Mat<zz_p>& Zp_out, const Mat<zz_p>& Yp_in, const Mat<zz_p>& Zp_in,
			     const zz_p& u1, const zz_p& v1, const zz_p& rho){


  Yp_out = Yp_in;
  Zp_out = Zp_in;

  long alpha = Yp_in.NumCols();
  long m = Yp_in.NumRows();
  long n = Zp_in.NumRows();
  long N = min(m, n);
  long p = zz_p::modulus();
  long step = alpha;

  Vec<zz_p> inverses_rho;
  Vec<zz_p> inverses_u1_v1;
  Vec<zz_p> inverses_u1_u1;
  Vec<zz_p> inverses_v1_v1;
  Vec<mulmod_precon_t> inverses_rho_pre;
  Vec<mulmod_precon_t> inverses_u1_v1_pre;
  Vec<mulmod_precon_t> inverses_u1_u1_pre;
  Vec<mulmod_precon_t> inverses_v1_v1_pre;

  inverse_powers(inverses_rho, rho, m);
  prepare_inverses_cauchy(inverses_u1_v1, u1, v1, rho, m, n); // 1/(u1-v1 rho^(-m+1)) ... 1/(u1-v1 rho^(n-1))
  // so i_u1_v1[i] = 1/1/(u1-v1 rho^(-m+1+i))   0 <= i < m+n-1
  // so i_u1_v1[i-j+m-1] = 1/1/(u1-v1 rho^(i-j))   0 <= i-j+m-1 < m+n-1  <=> i-n < j <= i+m-1
  prepare_inverses_cauchy(inverses_u1_u1, u1, u1, rho, N, 0); // 1/(u1-u1 rho^(-N+1)) ... 1/(u1-u1 rho^(-1))
  // so i_u1_u1[i] = 1/(u1-u1 rho^(-N+1+i))    0 <= i < N-1
  // so i_u1_u1[i-j+N-1]= 1/(u1-u1 rho^(i-j))   0 <= i-j+N-1 < N-1   <=>  i < j <= i+N-1
  prepare_inverses_cauchy(inverses_v1_v1, v1, v1, rho, N, 0); // 1/(v1-v1 rho^(-N+1)) ... 1/(v1-v1 rho^(-1))
  precomp(inverses_rho_pre, inverses_rho);
  precomp(inverses_u1_v1_pre, inverses_u1_v1);
  precomp(inverses_u1_u1_pre, inverses_u1_u1);
  precomp(inverses_v1_v1_pre, inverses_v1_v1);

  Mat<zz_p> gk, hk, tgk, thk, A, B, d, invd, neg_t_invd, mulG, mulH;

  bool early_exit = false;

  long k;
  for (k = 0; k < N; k += step){

    step = min(step, N-k);
    gk.SetDims(step, alpha);
    hk.SetDims(step, alpha);
    d.SetDims(step, step);
    neg_t_invd.SetDims(step, step);

    for (long i = 0; i < step; i++){
      gk[i] = Yp_out[i + k];
      hk[i] = Zp_out[i + k];
    }
    tgk = transpose(gk);
    thk = transpose(hk);

    mul(A, Yp_out, thk);
    mul(B, Zp_out, tgk);

    for (long i = 0; i < k; i++)
      for (long j = 0; j < step; j++){
    	// // (y[i]-y[k+j]);
	long ell = MulModPrecon((-A[i][j])._zz_p__rep, inverses_rho[k+j]._zz_p__rep, p, inverses_rho_pre[k+j]);
	A[i][j]._zz_p__rep = MulModPrecon(ell, inverses_v1_v1[i-(k+j)+N-1]._zz_p__rep, p, inverses_v1_v1_pre[i-(k+j)+N-1]);

    	// // (x[i]-x[k+j])
	long u = MulModPrecon((-B[i][j])._zz_p__rep, inverses_rho[k+j]._zz_p__rep, p, inverses_rho_pre[k+j]);
	B[i][j]._zz_p__rep = MulModPrecon(u, inverses_u1_u1[i-(k+j)+N-1]._zz_p__rep, p, inverses_u1_u1_pre[i-(k+j)+N-1]);  
      }

    for (long i = k; i < N; i++)
      for (long j = 0; j < step; j++){
    	// // (x[i]-y[k+j])
	long ell = MulModPrecon(A[i][j]._zz_p__rep, inverses_rho[i]._zz_p__rep, p, inverses_rho_pre[i]);
	A[i][j]._zz_p__rep = MulModPrecon(ell, inverses_u1_v1[k+j-i+m-1]._zz_p__rep, p, inverses_u1_v1_pre[k+j-i+m-1]);
    	// // -(x[k+j]-y[i]);
	long u = MulModPrecon((-B[i][j])._zz_p__rep, inverses_rho[k+j]._zz_p__rep, p, inverses_rho_pre[k+j]);
	B[i][j]._zz_p__rep = MulModPrecon(u, inverses_u1_v1[-(k+j)+i+m-1]._zz_p__rep, p, inverses_u1_v1_pre[-(k+j)+i+m-1]);
      }

    for (long i = 0; i < step; i++)
      for (long j = 0; j < step; j++)
	d[i][j] = -(A[i + k][j]);

    Mat<zz_p> d_cp = d;
    long rk = gauss(d_cp);

    if (rk < step){
      d_cp.SetDims(rk, rk);
      for (long i = 0; i < rk; i++)
	for (long j = 0; j < rk; j++)
	  d_cp[i][j] = d[i][j];
      long rk2 = gauss(d_cp);
      if (rk2 != rk){
	cerr << "the square block does not have grp\n";
	return -1;
      }
      k -= rk;  // go back
      step = rk;
      early_exit = true;
      continue;
    }

    inv(invd, d);

    for (long i = 0; i < step; i++)
      for (long j = 0; j < step; j++)
	neg_t_invd[i][j] = -invd[j][i];

    mul(mulG, invd, gk);
    mul(mulH, neg_t_invd, hk);
    mul(A, A, mulG);
    mul(B, B, mulH);

    add(Yp_out, Yp_out, A);
    add(Zp_out, Zp_out, B);

    for (long i = k; i < k + step; i++){
      Yp_out[i] = mulG[i-k];
      Zp_out[i] = mulH[i-k];
    }
    
    if (early_exit == true){
      Mat<zz_p> Y2, Z2;
      Y2.SetDims(m-(k+step), alpha);
      Z2.SetDims(n-(k+step), alpha);
      for (long i = k+step; i < m; i++)
	Y2[i-(k+step)] = Yp_out[i];
      for (long i = k+step; i < n; i++)
	Z2[i-(k+step)] = Zp_out[i];
      if (IsZero(Y2*transpose(Z2)))
	return -1;
      
      Yp_out.SetDims(k+step, alpha);
      Zp_out.SetDims(k+step, alpha);
      return k+step;
    }
  }
  
  Yp_out.SetDims(k, alpha);
  Zp_out.SetDims(k, alpha);
  return k;
}


/*------------------------------------------------------------*/
/* inverts a special-Cauchy-like matrix                       */
/* assumes generic rank profile                               */
/* the generators are Yp, Zp                                  */
/* alpha = displacement rank                                  */
/* overwrites Yp, Zp                                          */
/* i0, j0 are the initial indices in i, j resp.               */
/* i1, j1 are the end indices in i, j resp.                   */  
/*------------------------------------------------------------*/
static long invert_raw(Mat<zz_p>& Yp_out, Mat<zz_p>& Zp_out, const Mat<zz_p>& Yp_in, const Mat<zz_p>& Zp_in, 
		       const zz_p& u1, const zz_p& v1, const zz_p& rho){

  long alpha = Yp_in.NumCols();
  long m = Yp_in.NumRows();
  long n = Zp_in.NumRows();
  long N = min(m, n);

  long *Yp = new long[alpha*m];
  long *Zp = new long[alpha*n];

  long idx;

  idx = 0;
  for (long i = 0; i < m; i++){
    const zz_p *row = Yp_in[i].elts();
    for (long j = 0; j < alpha; j++)
      Yp[idx++] = row[j]._zz_p__rep;
  }
  idx = 0;
  for (long i = 0; i < n; i++){
    const zz_p *row = Zp_in[i].elts();
    for (long j = 0; j < alpha; j++)
      Zp[idx++] = row[j]._zz_p__rep;
  }


  const long p = zz_p::modulus();
  mulmod_precon_t *Ypre = new mulmod_precon_t[alpha];
  mulmod_precon_t *Zpre = new mulmod_precon_t[alpha];


  Vec<zz_p> inverses_rho;
  Vec<zz_p> inverses_u1_v1;
  Vec<zz_p> inverses_u1_u1;
  Vec<zz_p> inverses_v1_v1;
  Vec<mulmod_precon_t> inverses_rho_pre;
  Vec<mulmod_precon_t> inverses_u1_v1_pre;
  Vec<mulmod_precon_t> inverses_u1_u1_pre;
  Vec<mulmod_precon_t> inverses_v1_v1_pre;

  inverse_powers(inverses_rho, rho, m);
  prepare_inverses_cauchy(inverses_u1_v1, u1, v1, rho, m, n); // 1/(u1-v1 rho^(-m+1)) ... 1/(u1-v1 rho^(n-1))
  // so i_u1_v1[i] = 1/1/(u1-v1 rho^(-m+1+i))   0 <= i < m+n-1
  // so i_u1_v1[i-j+m-1] = 1/1/(u1-v1 rho^(i-j))   0 <= i-j+m-1 < m+n-1  <=> i-n < j <= i+m-1
  prepare_inverses_cauchy(inverses_u1_u1, u1, u1, rho, N, 0); // 1/(u1-u1 rho^(-N+1)) ... 1/(u1-u1 rho^(-1))
  // so i_u1_u1[i] = 1/(u1-u1 rho^(-N+1+i))    0 <= i < N-1
  // so i_u1_u1[i-j+N-1]= 1/(u1-u1 rho^(i-j))   0 <= i-j+N-1 < N-1   <=>  i < j <= i+N-1
  prepare_inverses_cauchy(inverses_v1_v1, v1, v1, rho, N, 0); // 1/(v1-v1 rho^(-N+1)) ... 1/(u1-u1 rho^(-1))
  precomp(inverses_rho_pre, inverses_rho);
  precomp(inverses_u1_v1_pre, inverses_u1_v1);
  precomp(inverses_u1_u1_pre, inverses_u1_u1);
  precomp(inverses_v1_v1_pre, inverses_v1_v1);

  const mulmod_t pinv = zz_p::ModulusInverse();
  long *tmpYk = Yp;
  long *tmpZk = Zp;

  for (long k = 0; k < N; k++){

    long *tmpYi = Yp;
    long *tmpZi = Zp;

    for (long j = 0; j < alpha; j++){
      Ypre[j] = PrepMulModPrecon(tmpYk[j], p, pinv);
      Zpre[j] = PrepMulModPrecon(tmpZk[j], p, pinv);
    }

    // computes the top-left entry in the current submatrix
    long rd = 0;
    for (long j = 0; j < alpha; j++)
      rd = AddMod(rd, MulMod(tmpZk[j], tmpYk[j], p, pinv), p);

    // early exit if rank(M) < n
    if (rd == 0){

      long *tmpY = tmpYk;
      for (long a = k; a < m; a++){
    	long *tmpZ = tmpZk;
    	for (long b = k; b < n; b++){
    	  long tmp = 0;
    	  for (long j = 0; j < alpha; j++)
    	    tmp = AddMod(tmp, MulMod(tmpZ[j], tmpY[j], p, pinv), p);
    	  if (tmp != 0) {
    	    delete[] Ypre;
    	    delete[] Zpre;
	    delete[] Yp;
	    delete[] Zp;
    	    return -1;  // not generic rank profile
    	  }
    	  tmpZ += alpha;
    	}
    	tmpY += alpha;
      }
      
      Yp_out.SetDims(k, alpha);
      idx = 0;
      for (long i = 0; i < k; i++){
	zz_p *row = Yp_out[i].elts();
	for (long j = 0; j < alpha; j++)
	  row[j] = Yp[idx++];
      }
      Zp_out.SetDims(k, alpha);
      idx = 0;
      for (long i = 0; i < k; i++){
	zz_p *row = Zp_out[i].elts();
	for (long j = 0; j < alpha; j++)
	  row[j] = Zp[idx++];
      }
      delete[] Yp;
      delete[] Zp;
      delete[] Ypre;
      delete[] Zpre;
      return k; // OK
    }

    rd = MulModPrecon(rd, inverses_rho[k]._zz_p__rep, p, inverses_rho_pre[k]);
    rd = MulModPrecon(rd, inverses_u1_v1[m-1]._zz_p__rep, p, inverses_u1_v1_pre[m-1]);

    zz_p d(rd, INIT_LOOP_HOLE);
    zz_p id = 1/d;
    zz_p mid = -id;

    const long rid = rep(id);
    const long rmid = rep(mid);
    mulmod_precon_t pd = PrepMulModPrecon(rd, p, pinv);
    mulmod_precon_t prmid = PrepMulModPrecon(rmid, p, pinv);

    for (long i = 0; i < k; i++){

      long ell = 0, u = 0;
      for (long j = 0; j < alpha; j++){
      	ell = AddMod(ell, MulModPrecon(tmpYi[j], tmpZk[j], p, Zpre[j]), p);
      	u = AddMod(u, MulModPrecon(tmpZi[j], tmpYk[j], p, Ypre[j]), p);
      }

      ell = MulModPrecon(ell, rmid, p, prmid);
      ell = MulModPrecon(ell, inverses_rho[k]._zz_p__rep, p, inverses_rho_pre[k]);
      ell = MulModPrecon(ell, inverses_v1_v1[i-k+N-1]._zz_p__rep, p, inverses_v1_v1_pre[i-k+N-1]);  // 1/(yk-yi)

      u = MulModPrecon(u, rmid, p, prmid);
      u = MulModPrecon(u, inverses_rho[k]._zz_p__rep, p, inverses_rho_pre[k]);
      u = MulModPrecon(u, inverses_u1_u1[i-k+N-1]._zz_p__rep, p, inverses_u1_u1_pre[i-k+N-1]);  // 1/(xk-xi)

      for (long j = 0; j < alpha; j++){
      	tmpYi[j] = SubMod(tmpYi[j], MulModPrecon(ell, tmpYk[j], p, Ypre[j]), p);
      	tmpZi[j] = AddMod(tmpZi[j], MulModPrecon(u, tmpZk[j], p, Zpre[j]), p);
      }

      tmpYi += alpha;
      tmpZi += alpha;
    }

    for (long j = 0; j < alpha; j++){
      tmpYk[j] = MulModPrecon(rmid, tmpYk[j], p, Ypre[j]);
      Ypre[j] = PrepMulModPrecon(tmpYk[j], p, pinv);
      tmpZk[j] = MulModPrecon(rid, tmpZk[j], p, Zpre[j]);
      Zpre[j] = PrepMulModPrecon(tmpZk[j], p, pinv);
    }

    tmpYi += alpha;
    tmpZi += alpha;

    for (long i = k+1; i < m; i++){

      long ell = 0;
      for (long j = 0; j < alpha; j++)
	ell = AddMod(ell, MulModPrecon(tmpYi[j], tmpZk[j], p, Zpre[j]), p);

      ell = MulModPrecon(ell, rd, p, pd);  
      ell = MulModPrecon(ell, inverses_rho[i]._zz_p__rep, p, inverses_rho_pre[i]);
      ell = MulModPrecon(ell, inverses_u1_v1[k-i+m-1]._zz_p__rep, p, inverses_u1_v1_pre[k-i+m-1]);

      for (long j = 0; j < alpha; j++)
      	tmpYi[j] = AddMod(tmpYi[j], MulModPrecon(ell, tmpYk[j], p, Ypre[j]), p);

      tmpYi += alpha;
    }

    for (long i = k+1; i < n; i++){

      long u = 0;
      for (long j = 0; j < alpha; j++)
        u = AddMod(u, MulModPrecon(tmpZi[j], tmpYk[j], p, Ypre[j]), p);
      
      // 1/(xk-yi)
      u = MulModPrecon(u, rd, p, pd); 
      u = MulModPrecon(u, inverses_rho[k]._zz_p__rep, p, inverses_rho_pre[k]);
      u = MulModPrecon(u, inverses_u1_v1[i-k+m-1]._zz_p__rep, p, inverses_u1_v1_pre[i-k+m-1]);

      for (long j = 0; j < alpha; j++)
      	tmpZi[j] = AddMod(tmpZi[j], MulModPrecon(u, tmpZk[j], p, Zpre[j]), p);

      tmpZi += alpha;
    }


    tmpYk += alpha;
    tmpZk += alpha;
  }

  Yp_out.SetDims(N, alpha);
  idx = 0;
  for (long i = 0; i < N; i++){
    zz_p *row = Yp_out[i].elts();
    for (long j = 0; j < alpha; j++)
      row[j] = Yp[idx++];
  }
  Zp_out.SetDims(N, alpha);
  idx = 0;
  for (long i = 0; i < N; i++){
    zz_p *row = Zp_out[i].elts();
    for (long j = 0; j < alpha; j++)
      row[j] = Zp[idx++];
  }
  delete[] Yp;
  delete[] Zp;
  delete[] Ypre;
  delete[] Zpre;

  return N;
}

/*---------------------------------------------------*/
/* iCL = CL^(-1), represented for the dual operator  */
/* O(alpha^2 M(n)) algorithm                         */
/*---------------------------------------------------*/
long invert_rec(Mat<zz_p>& Yp_out, Mat<zz_p>& Zp_out, 
		     const Mat<zz_p>& Yp_in, const Mat<zz_p>& Zp_in, 
		     const zz_p& u1, const zz_p& v1, const zz_p& rho,
		     const long thresh){

  long m = Yp_in.NumRows();
  long n = Zp_in.NumRows();
  long alpha = Yp_in.NumCols();
  long N = min(m, n);

  if (N < thresh){
    if (alpha < 15)
      return invert_raw(Yp_out, Zp_out, Yp_in, Zp_in, u1, v1, rho);
    else
      return invert_block_raw(Yp_out, Zp_out, Yp_in, Zp_in, u1, v1, rho);
  }

  long N1 = (N+1)/2;

  Mat<zz_p> Y1, Z1;
  Y1.SetDims(N1, alpha);
  Z1.SetDims(N1, alpha);
  for (long i = 0; i < N1; i++){
    Y1[i] = Yp_in[i];
    Z1[i] = Zp_in[i];
  }

  // Yp_out and Zp_out have size (r x alpha)
  long r = invert_rec(Yp_out, Zp_out, Y1, Z1, u1, v1, rho, thresh);

  // case 1: we found that we do not have grp. abort.
  if (r == -1)
    return -1;

  Y1.SetDims(r, alpha);
  Z1.SetDims(r, alpha);

  Mat<zz_p> G2;
  G2.SetDims(m - r, alpha);
  for (long i = 0; i < m - r; i++)
    G2[i] = Yp_in[i + r];

  Mat<zz_p> H2;
  H2.SetDims(n - r, alpha);
  for (long i = 0; i < n - r; i++)
    H2[i] = Zp_in[i + r];

  Mat<zz_p> tmp, tmp2;
  lzz_p_cauchy_like_geometric(G2, Z1, u1*power(rho, r), v1, rho).mul_right(tmp, Yp_out);
  // mul(tmp2, G2, Z1, Yp_out, u1*power(rho, r), v1, rho);
  // assert (tmp == tmp2);
  G2 += tmp;
  lzz_p_cauchy_like_geometric(H2, Y1, v1*power(rho, r), u1, rho).mul_right(tmp, Zp_out);
  // mul(tmp2, H2, Y1, Zp_out, v1*power(rho, r), u1, rho); 
  // assert (tmp == tmp2);
  H2 += tmp;

  // case 2: we found the rank of the top-left, which was grp. check all the rest.
  if (r < N1){
    // check if G2.H2^t = 0 by computing random_vec.G2.H2^t*random_vec
    Vec<zz_p> check;
    random(check, n-r);
    check = check*H2;
    check = G2*check;
    for (long i = 0; i < m-r; i++)
      if (check[i] != 0)
	return -1;
    return r;

  }

  // case 3: the top block was invertible.
  // now r=N1
  // iG2 and iH2 have size (s x alpha)
  Mat<zz_p> iG2, iH2;
  long s = invert_rec(iG2, iH2, G2, H2, u1*power(rho, N1), v1*power(rho, N1), rho, thresh);

  if (s == -1)
    return -1;

  G2.SetDims(s, alpha);
  H2.SetDims(s, alpha);
  
  lzz_p_cauchy_like_geometric(Yp_out, H2, v1, v1*power(rho, N1), rho).mul_right(tmp, iG2);
  // mul(tmp2, Yp_out, H2, iG2, v1, v1*power(rho, N1), rho);
  // assert (tmp2 == tmp);
  Yp_out += tmp;
  lzz_p_cauchy_like_geometric(Zp_out, G2, u1, u1*power(rho, N1), rho).mul_right(tmp, iH2);
  // mul(tmp2, Zp_out, G2, iH2, u1, u1*power(rho, N1), rho);
  // assert (tmp2 == tmp);
  Zp_out += tmp;
  
  Yp_out.SetDims(N1 + s, alpha);
  Zp_out.SetDims(N1 + s, alpha);

  for (long i = 0; i < s; i++){
    Yp_out[N1 + i] = iG2[i];
    Zp_out[N1 + i] = iH2[i];
  }

  return N1+s;
}


/*---------------------------------------------------*/
/* iCL = CL^(-1), represented for the dual operator  */
/* O(alpha n^2) algorithm                            */
/*---------------------------------------------------*/
long invert_direct(lzz_p_cauchy_like_geometric& Cinv,
		   const lzz_p_cauchy_like_geometric& CL){

  Mat<zz_p> Yp_out, Zp_out;

  long r = invert_raw(Yp_out, Zp_out, CL.G, CL.H, CL.C.u1, CL.C.v1, CL.C.rho);

  if (r == -1)
    return -1;

  Cinv = lzz_p_cauchy_like_geometric(Yp_out, Zp_out, CL.C.v1, CL.C.u1, CL.C.rho);
  return r;
}

/*---------------------------------------------------*/
/* iCL = CL^(-1), represented for the dual operator  */
/* O(alpha^{omega-2} n^2) algorithm                  */
/*---------------------------------------------------*/
long invert_block(lzz_p_cauchy_like_geometric& Cinv,
		      const lzz_p_cauchy_like_geometric& CL){

  Mat<zz_p> Yp_out, Zp_out;

  long r;
  if (CL.NumGens() > min(CL.NumRows(), CL.NumCols()))
    r = invert_raw(Yp_out, Zp_out, CL.G, CL.H, CL.C.u1, CL.C.v1, CL.C.rho);
  else
    r = invert_block_raw(Yp_out, Zp_out, CL.G, CL.H, CL.C.u1, CL.C.v1, CL.C.rho);

  if (r == -1)
    return -1;

  Cinv = lzz_p_cauchy_like_geometric(Yp_out, Zp_out, CL.C.v1, CL.C.u1, CL.C.rho);

  return r;
}

/*---------------------------------------------------*/
/* iCL = CL^(-1), represented for the dual operator  */
/* quadratic in n                                    */
/*---------------------------------------------------*/
long invert(lzz_p_cauchy_like_geometric& Cinv,
	    const lzz_p_cauchy_like_geometric& CL){

  if (CL.NumGens() < 15)
    return invert_direct(Cinv, CL);
  else
    return invert_block(Cinv, CL);
}

/*---------------------------------------------------*/
/* iCL = CL^(-1), represented for the dual operator  */
/* O(alpha^2 M(n)log(n)) algorithm                   */
/*---------------------------------------------------*/
long invert_fast(lzz_p_cauchy_like_geometric& Cinv,
		 const lzz_p_cauchy_like_geometric& CL, const long thresh){

  long do_thresh;
  if (thresh == -1)
    do_thresh = threshold(CL.NumGens());
  else
    do_thresh = thresh;

  Mat<zz_p> Yp_out, Zp_out;

  long r = invert_rec(Yp_out, Zp_out, CL.G, CL.H, CL.C.u1, CL.C.v1, CL.C.rho, do_thresh);

  if (r == -1)
    return -1;

  Cinv = lzz_p_cauchy_like_geometric(Yp_out, Zp_out, CL.C.v1, CL.C.u1, CL.C.rho);
  return r;
}

