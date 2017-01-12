#include <NTL/lzz_pX.h>
#include <NTL/mat_lzz_p.h>

#include "lzz_p_extra.h"
#include "lzz_pX_CRT.h"
#include "ZZ_pX_CRT.h"
#include "lzz_pX_mosaic_hankel.h"

NTL_CLIENT

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/* Mosaic Hankel matrices:                            */
/* block matrix where each block is Hankel            */
/*----------------------------------------------------*/
/*----------------------------------------------------*/

/*----------------------------------------------------*/
/* turns M into a dense matrix                        */
/*----------------------------------------------------*/
void to_dense(Mat<zz_p>& Mdense, const mosaic_hankel& M){
  long n = M.NumRows();
  long m = M.NumCols();
  Mdense.SetDims(n, m);

  for (long i = 0; i < n; i++)
    for (long j = 0; j < m; j++)
      Mdense[i][j] = M(i,j);
}

/*----------------------------------------------------*/
/* right multiplication                               */
/*----------------------------------------------------*/
void mul_right(Vec<zz_p>& res, const mosaic_hankel& M, const Vec<zz_p>& input){
  long n = M.NumRows();
  res.SetLength(n);
  for (long i = 0; i < n; i++)
    res[i] = 0;

  long nb = M.NumBlockRows();
  long mb = M.NumBlockCols();

  long maxn = 0;
  for (long i = 0; i < nb; i++)
    maxn = max(maxn, M.data[i][0].NumRows());

  long maxm = 0;
  for (long i = 0; i < mb; i++)
    maxm = max(maxm, M.data[0][i].NumCols());
  
  Vec<zz_p> tmp_in, tmp_out;
  tmp_in.SetLength(maxm);
  tmp_out.SetLength(maxn);

  long jdx = 0;
  for (long j = 0; j < mb; j++){
    for (long ell = 0; ell < M.data[0][j].NumCols(); ell++)
      tmp_in[ell] = input[jdx + ell];

    long idx = 0;
    for (long i = 0; i < nb; i++){
      mul_right(tmp_out, M.data[i][j], tmp_in);
      for (long ell = 0; ell < M.data[i][0].NumRows(); ell++)
	res[idx + ell] += tmp_out[ell];
      idx += M.data[i][0].NumRows();
    }
    jdx += M.data[0][j].NumCols();
  }
}

/*----------------------------------------------------*/
/* left multiplication                                */
/*----------------------------------------------------*/
void mul_left(Vec<zz_p>& res, const mosaic_hankel& M, const Vec<zz_p>& input){
  long m = M.NumCols();
  res.SetLength(m);
  for (long i = 0; i < m; i++)
    res[i] = 0;

  long nb = M.NumBlockRows();
  long mb = M.NumBlockCols();

  long maxn = 0;
  for (long i = 0; i < nb; i++)
    maxn = max(maxn, M.data[i][0].NumRows());

  long maxm = 0;
  for (long i = 0; i < mb; i++)
    maxm = max(maxm, M.data[0][i].NumCols());

  Vec<zz_p> tmp_in, tmp_out;
  tmp_in.SetLength(maxn);
  tmp_out.SetLength(maxm);

  long idx = 0;
  for (long i = 0; i < nb; i++){
    for (long ell = 0; ell < M.data[i][0].NumRows(); ell++)
      tmp_in[ell] = input[idx + ell];

    long jdx = 0;
    for (long j = 0; j < mb; j++){
      mul_left(tmp_out, M.data[i][j], tmp_in);
      for (long ell = 0; ell < M.data[0][j].NumCols(); ell++)
	res[jdx + ell] += tmp_out[ell];
      jdx += M.data[0][j].NumCols();
    }
    idx += M.data[i][0].NumRows();
  }

}

/*----------------------------------------------------*/
/* access to particular rows and columns              */
/*----------------------------------------------------*/
void first_column_of_block(Vec<zz_p>& res, long i, const mosaic_hankel& M){
  long n = M.NumRows();
  long nb = M.NumBlockRows();

  res.SetLength(n);

  long ind = 0;
  for (long r = 0; r < nb; r++){
    const zz_p* dat = M.data[r][i].data_rev.elts();
    for (long j = 0; j < M.NumRows_of_block(r); j++)
      res[ind++] = dat[j];
  }
}

/*----------------------------------------------------*/
/* access to particular rows and columns              */
/*----------------------------------------------------*/
void last_column_of_block(Vec<zz_p>& res, long i, const mosaic_hankel& M){
  long n = M.NumRows();
  long nb = M.NumBlockRows();

  res.SetLength(n);

  long ind = 0;
  for (long r = 0; r < nb; r++){
    const zz_p* dat = M.data[r][i].data_rev.elts();
    long shift = M.NumCols_of_block(i)-1;
    for (long j = 0; j < M.NumRows_of_block(r); j++)
      res[ind++] = dat[j + shift];
  }
}

/*----------------------------------------------------*/
/* access to particular rows and columns              */
/*----------------------------------------------------*/
void first_row_of_block(Vec<zz_p>& res, long i, const mosaic_hankel& M){
  long m = M.NumCols();
  long mb = M.NumBlockCols();

  res.SetLength(m);

  long ind = 0;
  for (long r = 0; r < mb; r++){
    const zz_p* dat = M.data[i][r].data_rev.elts();
    for (long j = 0; j < M.NumCols_of_block(r); j++)
      res[ind++] = dat[j];
  }
}

/*----------------------------------------------------*/
/* access to particular rows and columns              */
/*----------------------------------------------------*/
void last_row_of_block(Vec<zz_p>& res, long i, const mosaic_hankel& M){
  long m = M.NumCols();
  long mb = M.NumBlockCols();

  res.SetLength(m);

  long ind = 0;
  for (long r = 0; r < mb; r++){
    const zz_p* dat = M.data[i][r].data_rev.elts();
    long shift = M.NumRows_of_block(i)-1;
    for (long j = 0; j < M.NumCols_of_block(r); j++)
      res[ind++] = dat[j+shift];
  }
}

/*----------------------------------------------------*/
/* G, H such that Z1 M - Z0^t M = G H^t               */
/*----------------------------------------------------*/
void generators(Mat<zz_p>& G, Mat<zz_p>& H, const mosaic_hankel& M){
  long n = M.NumRows();
  long m = M.NumCols();
  long alpha_r = M.NumBlockRows();
  long alpha_c = M.NumBlockCols();
  long alpha = alpha_r + alpha_c;

  G.SetDims(n, alpha);
  H.SetDims(m, alpha);

  Vec<zz_p> tmp_row, tmp_col;
  tmp_row.SetLength(m);
  tmp_col.SetLength(n);

  long idx;

  idx = 0;
  for (long i = 0; i < alpha_c; i++){
    first_column_of_block(tmp_col, i, M);
    G[0][i] = tmp_col[n-1];
    for (long j = 1; j < n; j++)
      G[j][i] = tmp_col[j-1];

    if (i > 0){
      last_column_of_block(tmp_col, i-1, M);
      for (long j = 0; j < n; j++)
    	G[j][i] = G[j][i] - tmp_col[j];
    }
  }
  for (long i = 0; i < alpha_r; i++){
    for (long j = 0; j < n; j++)
      G[j][i+alpha_c] = 0;
    G[idx][i+alpha_c] = 1;
    idx += M.NumRows_of_block(i);
  }


  idx = 0;
  for (long i = 0; i < alpha_c; i++){
    for (long j = 0; j < m; j++)
      H[j][i] = 0;
    H[idx][i] = 1;
    idx += M.NumCols_of_block(i);
  }
  for (long i = 0; i < alpha_r; i++){
    if (i == 0)
      last_row_of_block(tmp_row, alpha_r-1, M);
    else
      last_row_of_block(tmp_row, i-1, M);
    for (long j = 1; j < m; j++)
      H[j][i+alpha_c] = tmp_row[j];
    first_row_of_block(tmp_row, i, M);
    for (long j = 1; j < m; j++)
      H[j][i+alpha_c] = H[j][i+alpha_c] - tmp_row[j-1];
    long jdx = 0;
    for (long j = 0; j < alpha_c; j++){
      H[jdx][i+alpha_c] = 0;
      jdx += M.NumCols_of_block(j);
    }
  }
}

/*------------------------------------------------------------------*/
/* finds c such that                                                */
/* - c != 0                                                         */
/* - a^i - c a^j != 0 for 0 <= i < n and 0 <= j < m                 */
/*------------------------------------------------------------------*/
static 
void find_c(zz_p& c, const zz_p& a, long n, long m){
  Vec<zz_p> pow_a, pow_inva;
  pow_a.SetLength(n);
  pow_inva.SetLength(m);

  pow_a[0] = to_zz_p(1);
  pow_inva[0] = to_zz_p(1);

  zz_p inva = 1/a;
  for (long i = 1; i < n; i++)
    pow_a[i] = a*pow_a[i-1];
  for (long i = 1; i < m; i++)
    pow_inva[i] = inva*pow_inva[i-1];

  bool done;
  do{
    c = random_zz_p();
    done = true;
    if (c == 0)
      done = false;
    for (long i = 0; i < n; i++)
      if (c == pow_a[i])
	done = false;
    for (long i = 0; i < m; i++)
      if (c == pow_inva[i])
	done = false;
  } while (done != true);
}

/*------------------------------------------------------------------*/
/* preconditions M                                                  */
/* builds the matrix CL = (D_e X_int) M (D_f Y_int)^t, where:       */
/* - X_int, Y_int are geometric interpolation                       */
/* - D_e, D_f are diagonal matrix built on vectors e and f          */
/* - CL is cauchy-like                                              */
/* - CL is expected to have generic rank profile                    */
/*                                                                  */
/* - X_int is built on (1,w,w^2,..)                                 */
/* - Y_int is built on (c,cw,cw^2,..)                               */
/* If we are over an FFT prime, w = root of unity                   */
/*------------------------------------------------------------------*/
void to_cauchy_grp(lzz_p_cauchy_like_geometric& CL, 
		   zz_pX_Multipoint_Geometric& X_int, zz_pX_Multipoint_Geometric& Y_int,
		   Vec<zz_p> &e, Vec<zz_p> &f,
		   const mosaic_hankel& M){


  Mat<zz_p> X, Y;
  Mat<zz_p> G, H;
  generators(G, H, M);
  lzz_p_cauchy_geometric C;
  zz_p a, b, c;
  long n = M.NumRows();
  long m = M.NumCols();

  if (max(m, n) > zz_p::modulus())
    LogicError("Field too small for preconditioning.");

  zz_pInfoT *info = zz_pInfo;
  if (info->p_info != NULL){  // FFT prime
    long k = NextPowerOfTwo(max(m, n));
    long order = 1L << k;
    long w = find_root_of_unity(zz_p::modulus(), 2*order);
    a = to_zz_p(w);
  }
  else
    element_of_order(a, max(m, n));

  b = a*a;
  find_c(c, b, n, m); 
  C = lzz_p_cauchy_geometric(to_zz_p(1), c, b, n, m);
  C.build_X_Y();

  X_int = C.X;
  Y_int = C.Y;

  long alpha = G.NumCols();
  X.SetDims(n, alpha+2);
  Y.SetDims(m, alpha+2);
      
  e.SetLength(n);
  for (long i = 0; i < n; i++)
    e[i] = random_zz_p();

  f.SetLength(m);
  for (long i = 0; i < m; i++)
    f[i] = random_zz_p();

  vec_zz_p tmp_v;
  for (long j = 0; j < alpha; j++){
    zz_pX tmp_p;
    tmp_p.rep.SetLength(n);
    zz_p* coef_p = tmp_p.rep.elts();
    for (long i = 0; i < n; i++)
      coef_p[i] = G[i][j];
    tmp_p.normalize();
    X_int.evaluate(tmp_v, tmp_p);
    for (long i = 0; i < n; i++)
      X[i][j] = tmp_v[i] * e[i];
  }

  zz_p tmp_z = to_zz_p(1);
  for (long i = 0; i < n; i++){
    X[i][alpha] = e[i]*(power(tmp_z, n)-1);
    tmp_z = tmp_z * b;
  }

  last_column_of_block(tmp_v, M.NumBlockCols()-1, M);
  zz_pX tmp_p;
  tmp_p.rep.SetLength(n);
  zz_p* coef_p = tmp_p.rep.elts();
  for (long i = 0; i < n; i++)
    coef_p[i] = tmp_v[i];
  tmp_p.normalize();
  X_int.evaluate(tmp_v, tmp_p);
  for (long i = 0; i < n; i++)
    X[i][alpha+1] = tmp_v[i] * e[i];

  vec_zz_p tmp_w;
  for (long j = 0; j < alpha; j++){
    zz_pX tmp_q;
    tmp_q.rep.SetLength(m);
    zz_p* coef_q = tmp_q.rep.elts();
    for (long i = 0; i < m; i++)
      coef_q[i] = H[i][j];
    tmp_q.normalize();
    Y_int.evaluate(tmp_w, tmp_q);
    for (long i = 0; i < m; i++)
      Y[i][j] = tmp_w[i] * f[i];
  }

  last_row_of_block(tmp_w, M.NumBlockRows()-1, M);
  zz_pX tmp_q;
  tmp_q.rep.SetLength(m);
  zz_p* coef_q = tmp_q.rep.elts();
  for (long i = 0; i < m; i++)
    coef_q[i] = tmp_w[i];
  tmp_q.normalize();
  Y_int.evaluate(tmp_w, tmp_q);
  for (long i = 0; i < m; i++)
    Y[i][alpha] = tmp_w[i] * f[i];

  tmp_z = c;
  for (long i = 0; i < m; i++){
    Y[i][alpha+1] = -f[i]*power(tmp_z, m);
    tmp_z = tmp_z * b;
  }

  CL = lzz_p_cauchy_like_geometric(X, Y, to_zz_p(1), c, b);
}
