#include <NTL/lzz_pX.h>
#include <NTL/mat_lzz_p.h>

#include "lzz_pX_CRT.h"
#include "lzz_p_toeplitz_like.h"

NTL_CLIENT

lzz_p_toeplitz_like::lzz_p_toeplitz_like(const Mat<zz_p>& G0, const Mat<zz_p>& H0){
  G = G0;
  H = H0;
  alpha = G.NumCols();
  m = G.NumRows();
  n = H.NumRows();
}

/*----------------------------------------------------*/
/* output = M^t * input                               */
/*----------------------------------------------------*/
void lzz_p_toeplitz_like::mul_left(Vec<zz_p>& output, const Vec<zz_p>& input){
  zz_pX f, res;
 
  for (long i = 0; i < m; i++)
    SetCoeff(f, i, input[i]);

  res = 0;

  for (long j = 0; j < alpha; j++){
    zz_pX t, g, h;
    for (long i = 0; i < m; i++)
      SetCoeff(g, i, G[m-1-i][j]);
    for (long i = 0; i < n; i++)
      SetCoeff(h, i, H[i][j]);
    t = trunc(g * f, m);
    res += (h * t);

  }

  output.SetLength(n);
  for (long i = 0; i < n; i++)
    output[i] = coeff(res, i + m - 1);
}

/*----------------------------------------------------*/
/* output = M * input                                 */
/*----------------------------------------------------*/
void lzz_p_toeplitz_like::mul_right(Vec<zz_p>& output, const Vec<zz_p>& input){
  zz_pX f, res;
 
  for (long i = 0; i < n; i++)
    SetCoeff(f, i, input[i]);

  res = 0;

  for (long j = 0; j < alpha; j++){
    zz_pX t, g, h;
    for (long i = 0; i < m; i++)
      SetCoeff(g, i, G[i][j]);
    for (long i = 0; i < n; i++)
      SetCoeff(h, i, H[n-1-i][j]);
    t = trunc(h * f, n); 
    res += (g * t);
  }

  output.SetLength(m);
  for (long i = 0; i < m; i++)
    output[i] = coeff(res, i + n - 1);
}


/*----------------------------------------------------*/
/*----------------------------------------------------*/
void lzz_p_toeplitz_like::mul_right_dac(Mat<zz_p>& output, const Mat<zz_p>& input){
  // TODO: check that we are FFT prime

  long beta = input.NumCols();

  output.SetDims(m, beta);
  for (long i = 0; i < m; i++)
    for (long j = 0; j < beta; j++)
      output[i][j] = 0;


  Vec<Mat<zz_p>> valA, valG, valH, valC;
  long h_col = 1;
  long nn = n;
  long shift = 0;

  // increments we use to extract polys from H and C
  Vec<long> deltas;

  while (h_col < alpha && nn > 4){

    long nn0 = nn / 2;
    long nn1 = nn - nn0;
    long eps = nn1 - nn0;
    long d = 3*nn1 - 2;
    long g_row = (m + nn1 - 1) / nn1;

    valA.SetLength(d);
    valG.SetLength(d);
    valH.SetLength(d);
    valC.SetLength(d);

    for (long i = 0; i < d; i++){
      valG[i].SetDims(g_row, alpha);
      valH[i].SetDims(alpha, h_col);
      valC[i].SetDims(h_col, beta);
    }

    zz_pX_Multipoint_CTFT ev(d);
    
    zz_pX tmp;
    Vec<zz_p> tmpval;
    tmpval.SetLength(d);

    for (long k = 0; k < alpha; k++){
      long idx;

      idx = 0;
      for (long j = 0; j < g_row; j++){
	tmp.rep.SetLength(nn1);
	long jj;
	for (jj = 0; jj < nn1 && idx < m; jj++)
	  tmp[jj] = G[idx++][k];
	for (; jj < nn1; jj++)
	  tmp[jj] = 0;
	tmp.normalize();
	ev.evaluate(tmpval, tmp);
	for (long jj = 0; jj < d; jj++)
	  valG[jj][j][k] = tmpval[jj];
      }

      idx = 0;
      for (long j = 0; j < h_col; j++){
	tmp.rep.SetLength(nn1);
	long jj;
	for (jj = 0; jj < nn1 && idx < n; jj++){
	  tmp[jj] = H[n - 1 - idx][k];
	  idx++;
	}
	for (; jj < nn1; jj++)
	  tmp[jj] = 0;
	tmp.normalize();
	ev.evaluate(tmpval, tmp);
	for (long jj = 0; jj < d; jj++)
	  valH[jj][k][j] = tmpval[jj];
	if (j < h_col - 1)
	  idx += nn0 + deltas[j];
      }
    }

    for (long k = 0; k < beta; k++){
      long idx;

      idx = 0;
      for (long j = 0; j < h_col; j++){
	tmp.rep.SetLength(nn1);
	long jj;
	for (jj = 0; jj < nn1 && idx < n; jj++)
	  tmp[jj] = input[idx++][k];
	for (; jj < nn1; jj++)
	  tmp[jj] = 0;
	tmp.normalize();
	ev.evaluate(tmpval, tmp);
	for (long jj = 0; jj < d; jj++)
	  valC[jj][h_col-j-1][k] = tmpval[jj];
	if (j < h_col - 1)
	  idx += nn0 + deltas[j];
      }
    }

    Mat<zz_p> tmpmat;
    for (long jj = 0; jj < d; jj++){
      mul(tmpmat, valG[jj], valH[jj]);
      mul(valA[jj], tmpmat, valC[jj]);
    }


    for (long k = 0; k < beta; k++){
      long idx = 0;
      for (long j = 0; j < g_row; j++){
	for (long jj = 0; jj < d; jj++)
	  tmpval[jj] = valA[jj][j][k];
	ev.interpolate(tmp, tmpval);

	long step = shift+idx-(n-1);
	long start = max(0, -step);
	long end = min(d, m-step);
	for (long jj = start; jj < end; jj++)
	  output[step+jj][k] += coeff(tmp, jj);
	idx += nn1;
      }
    }

    h_col *= 2;
    shift += nn1;
    nn = nn0;

    long len = deltas.length();
    Vec<long> cp = deltas;
    deltas.SetLength(2*len+1);
    deltas[0] = eps;
    for (long i = 0; i < len; i++){
      deltas[2*i+1] = cp[i];
      deltas[2*i+2] = eps;
    }
  }

  {
    long g_row = (m + nn - 1) / nn;

    Vec<Mat<zz_p>> valtmp;

    long d = 2*nn - 1;
    valA.SetLength(d); 
    valtmp.SetLength(d);
    valH.SetLength(d);
    valC.SetLength(d);

    for (long i = 0; i < d; i++){
      valG[i].SetDims(g_row, alpha);
      valH[i].SetDims(alpha, h_col);
      valC[i].SetDims(h_col, beta);
    }

    zz_pX_Multipoint_CTFT ev(d);
    
    zz_pX tmp;
    Vec<zz_p> tmpval;
    tmpval.SetLength(d);

    for (long k = 0; k < alpha; k++){
      long idx;

      idx = 0;
      for (long j = 0; j < g_row; j++){
	tmp.rep.SetLength(nn);
	long jj;
	for (jj = 0; jj < nn && idx < m; jj++)
	  tmp[jj] = G[idx++][k];
	for (; jj < nn; jj++)
	  tmp[jj] = 0;
	tmp.normalize();
	ev.evaluate(tmpval, tmp);
	for (long jj = 0; jj < d; jj++)
	  valG[jj][j][k] = tmpval[jj];
      }

      idx = 0;
      for (long j = 0; j < h_col; j++){
	tmp.rep.SetLength(nn);
	long jj;
	for (jj = 0; jj < nn && idx < n; jj++){
	  tmp[jj] = H[n - 1 - idx][k];
	  idx++;
	}
	for (; jj < nn; jj++)
	  tmp[jj] = 0;
	tmp.normalize();
	ev.evaluate(tmpval, tmp);
	for (long jj = 0; jj < d; jj++)
	  valH[jj][k][j] = tmpval[jj];
	if (j < h_col - 1)
	  idx += deltas[j];
      }
    }

    for (long k = 0; k < beta; k++){
      long idx;

      idx = 0;
      for (long j = 0; j < h_col; j++){
	tmp.rep.SetLength(nn);
	long jj;
	for (jj = 0; jj < nn && idx < n; jj++)
	  tmp[jj] = input[idx++][k];
	for (; jj < nn; jj++)
	  tmp[jj] = 0;
	tmp.normalize();
	ev.evaluate(tmpval, tmp);
	for (long jj = 0; jj < d; jj++)
	  valC[jj][h_col-j-1][k] = tmpval[jj];
	if (j < h_col - 1)
	  idx += deltas[j];
      }
    }

    for (long jj = 0; jj < d; jj++)
      mul(valtmp[jj], valH[jj], valC[jj]);

    for (long k = 0; k < beta; k++){
      for (long j = 0; j < alpha; j++){
	for (long jj = 0; jj < d; jj++)
	  tmpval[jj] = valtmp[jj][j][k];
	ev.interpolate(tmp, tmpval);
	ev.evaluate(tmpval, trunc(tmp, nn));
	for (long jj = 0; jj < d; jj++)
	  valtmp[jj][j][k] = tmpval[jj];
      }
    }

    for (long jj = 0; jj < d; jj++)
      mul(valA[jj], valG[jj], valtmp[jj]);

    for (long k = 0; k < beta; k++){
      long idx = 0;

      for (long j = 0; j < g_row; j++){
	for (long jj = 0; jj < d; jj++)
	  tmpval[jj] = valA[jj][j][k];

	ev.interpolate(tmp, tmpval);

	long step = shift+idx-(n-1);
	long start = max(0, -step);
	long end = min(d, m-step);
	for (long jj = start; jj < end; jj++)
	  output[step+jj][k] += coeff(tmp, jj);
	
	idx += nn;
      }
    }
  }

}


/*----------------------------------------------------*/
/* down-shift matrix of size n                        */
/*----------------------------------------------------*/
Mat<zz_p> do_Z(long n){
  Mat<zz_p> Z;
  Z.SetDims(n, n);
  for (long i = 0; i < n-1; i++)
    Z[i+1][i] = 1;
  return Z;
}

/*----------------------------------------------------*/
/* reconstructs M from its generators                 */
/*----------------------------------------------------*/
void lzz_p_toeplitz_like::to_dense(Mat<zz_p>& output){
  output.SetDims(m, n);
  Mat<zz_p> Zm = do_Z(m);
  Mat<zz_p> Zn = do_Z(n);
  Mat<zz_p> Gtmp = G;
  Mat<zz_p> Htmp = H;

  for (long i = 0; i < min(m, n); i++){
    output = output + Gtmp * transpose(Htmp);
    Gtmp = transpose(Zm) * Gtmp;
    Htmp = transpose(Zn) * Htmp;
  }
}

/*----------------------------------------------------*/
/* M - Z^t M Z = M - (M shifted up, left)             */
/*----------------------------------------------------*/
void phi_minus(Mat<zz_p>& output, const Mat<zz_p>& input){
  long m = input.NumRows();
  long n = input.NumCols();
  output.SetDims(m, n);
  for (long i = 0; i < m; i++)
    for (long j = 0; j < n; j++){
      output[i][j] = input[i][j];
      if (i < m-1 && j < n-1)
	output[i][j] -= input[i+1][j+1];
    }
}
