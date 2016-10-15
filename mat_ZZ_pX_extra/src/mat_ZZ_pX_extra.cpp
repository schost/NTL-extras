#include <NTL/matrix.h>
#include <NTL/ZZ_pX.h>
#include <NTL/lzz_pX.h>

#include "lzz_pX_CRT.h"
#include "mat_ZZ_pX_extra.h"
#include "mat_lzz_pX_extra.h"
#include "ZZ_pX_extra.h"
#include "ZZ_extra.h"
#include "magma_output.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/* magma output, without assign, using variable var           */
/*------------------------------------------------------------*/
void magma_output(const Mat<ZZ_pX>& a, const string & var){
  cout << "Matrix([";
  for (long i = 0; i < a.NumRows(); i++){
    cout << "[";
    for (long j = 0; j < a.NumCols(); j++){
      magma_output(a[i][j], var);
      if (j < a.NumCols()-1)
	cout << ", ";
    }
    cout << "]";
    if (i < a.NumRows()-1)
      cout << ", ";
  }
  cout << "])";
}

/*------------------------------------------------------------*/
/* magma assignment, using variable var, to name              */
/*------------------------------------------------------------*/
void magma_assign(const Mat<ZZ_pX>& a, const string & var, const string & name){
  cout << name << ":=";
  magma_output(a, var);
  cout << ";\n";
}

/*------------------------------------------------------------*/
/* random matrix of a given degree                            */
/*------------------------------------------------------------*/
void random_mat_ZZ_pX(Mat<ZZ_pX>& a, long n, long m, long d){
  a.SetDims(n, m);
  for (long i = 0; i < n; i++)
    for (long j = 0; j < m; j++)
      a[i][j] = random_ZZ_pX(d);
}

/*------------------------------------------------------------*/
/* maximum degree of the entries of a                         */
/*------------------------------------------------------------*/
long deg(const Mat<ZZ_pX>& a){
  long d = 0;
  long m = a.NumRows();
  long n = a.NumCols();

  for (long i = 0; i < m; i++)
    for (long j = 0; j < n; j++)
      d = max(deg(a[i][j]), d);
  return d;
}

/*------------------------------------------------------------*/
/* maximum size of the entries of a                           */
/*------------------------------------------------------------*/
long size(const Mat<ZZ_pX>& a){
  long d = 0;
  long m = a.NumRows();
  long n = a.NumCols();

  for (long i = 0; i < m; i++)
    for (long j = 0; j < n; j++)
      d = max(size(a[i][j]), d);
  return d;
}

/*------------------------------------------------------------*/
/*  direct multiplication                                     */
/*------------------------------------------------------------*/
void mul_direct(Mat<ZZ_pX>& X, const Mat<ZZ_pX>& A, const Mat<ZZ_pX>& B){  
   long n = A.NumRows();  
   long l = A.NumCols();  
   long m = B.NumCols();  
  
   if (l != B.NumRows())  
      Error("matrix mul: dimension mismatch");  
  
   X.SetDims(n, m);  
  
   long i, j, k;  
   ZZ_pX acc, tmp;  
  
   for (i = 1; i <= n; i++) {  
      for (j = 1; j <= m; j++) {  
         clear(acc);  
         for(k = 1; k <= l; k++) {  
	   mul(tmp, A(i,k), B(k,j));  
	   add(acc, acc, tmp);  
         }  
         X(i,j) = acc;  
      }  
   }  
}  

/*------------------------------------------------------------*/
/*  Waksman's algorithm                                       */
/*------------------------------------------------------------*/
void mul_waksman(Mat<ZZ_pX> &C, const Mat<ZZ_pX> &A, const Mat<ZZ_pX> &B){
  vec_ZZ_pX Arow, Acol, D, E;
  ZZ_pX val0, val1, val2, crow;

  long m = A.NumRows();
  long n = A.NumCols();
  long p = B.NumCols();

  C.SetDims(m, p);
  
  D.SetLength(p);
  E.SetLength(m);

  long np = n>>1;

  for (long j=1; j<=np; j++){
    const long j2=(j<<1)-1;
    
    for (long k=0; k<p; k++){
      C[0][k] += (A[0][j2-1]+B[j2][k]) * (A[0][j2]+B[j2-1][k]);
      D[k] += (A[0][j2-1]-B[j2][k]) * (A[0][j2]-B[j2-1][k]);
    }

    for (long l=1; l<m; l++){
      C[l][0] += (A[l][j2-1]+B[j2][0]) * (A[l][j2]+B[j2-1][0]);
      E[l] += (A[l][j2-1]-B[j2][0]) * (A[l][j2]-B[j2-1][0]);
    }

    for (long k=1; k<p; k++)
      for (long l=1; l<m; l++)
	C[l][k] += (A[l][j2-1]+B[j2][k]) * (A[l][j2]+B[j2-1][k]);
  }

  for (long l=1; l<m; l++){
    E[l] = (E[l]+C[l][0])/2;
    C[l][0] -= E[l];
  }
  val0 = (D[0]+C[0][0])/2;
  C[0][0] -= val0;
  for (long k=1; k<p; k++){
    val1 = (D[k]+C[0][k])/2;
    C[0][k] -= val1;
    val1 -= val0;
    for (long l=1; l<m; l++)
      C[l][k] -= val1 + E[l];
  }
  if ((n & 1) == 1){
    for (long l=0; l<m; l++)
      for (long k=0; k<p; k++)
  	C[l][k] += A[l][n-1]*B[n-1][k];
  }
}



void mul_CRT_CTFT(Mat<ZZ_pX>& C, const Mat<ZZ_pX>& A, const Mat<ZZ_pX>& B){
  // format a,b,c
  long a = A.NumRows();
  long b = A.NumCols();
  long c = B.NumCols();

  ZZ p = ZZ_p::modulus();
  Mat<ZZX> ZA, ZB, ZC;
  conv(ZA, A);
  conv(ZB, B);

  ZZ_pContext ctx;
  ctx.save();
  zz_pContext lctx;
  lctx.save();
  
  ZZ_p::init(p * (1 + _ntl_gsqrts(b)));
  Mat<ZZ_pX> A2, B2, C2;
  conv(A2, A);
  conv(B2, B);

  const ZZ_pFFTInfoT *FFTInfo = ZZ_p::GetFFTInfo();
  ZZ_pTmpSpaceT *TmpSpace = ZZ_p::GetTmpSpace();
  long nprimes = FFTInfo->NumPrimes;

  Vec<zz_pContext> vec_fft_contexts;
  vec_fft_contexts.SetLength(nprimes);
  for (long i = 0; i < nprimes; i++)
    vec_fft_contexts[i] = zz_pContext(INIT_FFT, i);

  long nA = deg(A2);
  long nB = deg(B2);
  long n = nA + nB + 1;
  Vec<zz_pX_Multipoint_CTFT> vec_CTFT;
  vec_CTFT.SetLength(nprimes);
  for (long i = 0; i < nprimes; i++){
    vec_fft_contexts[i].restore();
    vec_CTFT[i] = zz_pX_Multipoint_CTFT(n);
  }

  Vec<Vec<Mat<zz_p>>> valA, valB, valC;
  valA.SetLength(nprimes);
  valB.SetLength(nprimes);
  valC.SetLength(nprimes);

  Vec<Vec<long>> t;
  t.SetLength(nprimes);

  Vec<long> tmp;
  tmp.SetLength(nprimes);

  for (long i = 0; i < nprimes; i++){
    t[i].SetLength(n);
    valA[i].SetLength(n);
    valB[i].SetLength(n);
    valC[i].SetLength(n);
    for (long j = 0; j < n; j++){
      valA[i][j].SetDims(a, b);
      valB[i][j].SetDims(b, c);
    }
  }
  
  for (long u = 0; u < a; u++)
    for (long v = 0; v < b; v++){
      const ZZ_p *xx = A2[u][v].rep.elts();
      for (long j = 0; j <= deg(A2[u][v]); j++) {
	to_modular_rep(tmp, xx[j], FFTInfo, TmpSpace);
	for (long i = 0; i < nprimes; i++) 
	  t[i][j] = tmp[i];
      }

      for (long j = deg(A2[u][v])+1; j < n; j++) 
	for (long i = 0; i < nprimes; i++) 
	  t[i][j] = 0;

      for (long i = 0; i < nprimes; i++){
	vec_fft_contexts[i].restore();
	Vec<zz_p> val;
	vec_CTFT[i].evaluate(val, t[i].elts(), n);
	for (long j = 0; j < n; j++)
	  valA[i][j][u][v] = val[j];
      }
    }

  for (long u = 0; u < b; u++)
    for (long v = 0; v < c; v++){
      const ZZ_p *xx = B2[u][v].rep.elts();
      for (long j = 0; j <= deg(B2[u][v]); j++) {
	to_modular_rep(tmp, xx[j], FFTInfo, TmpSpace);
	for (long i = 0; i < nprimes; i++) 
	  t[i][j] = tmp[i];
      }

      for (long j = deg(B2[u][v])+1; j < n; j++) 
	for (long i = 0; i < nprimes; i++) 
	  t[i][j] = 0;

      for (long i = 0; i < nprimes; i++){
	vec_fft_contexts[i].restore();
	Vec<zz_p> val;
	vec_CTFT[i].evaluate(val, t[i].elts(), n);
	for (long j = 0; j < n; j++)
	  valB[i][j][u][v] = val[j];
      }
    }

  
  for (long i = 0; i < nprimes; i++){
    vec_fft_contexts[i].restore();
    for (long j = 0; j < n; j++)
      mul(valC[i][j], valA[i][j], valB[i][j]);
  }

  C2.SetDims(a, c);

  for (long u = 0; u < a; u++)
    for (long v = 0; v < c; v++){
      for(long i = 0; i < nprimes; i++){
	vec_fft_contexts[i].restore();
	Vec<zz_p> valf;
	valf.SetLength(n);
	zz_pX f;
	for (long j = 0; j < n; j++)
	  valf[j] = valC[i][j][u][v];
	vec_CTFT[i].interpolate(f, valf);
	long df = deg(f);
	zz_p * cf = f.rep.elts();
	for (long j = 0; j <= df; j++)
	  t[i][j] = cf[j]._zz_p__rep;
	for (long j = df+1; j < n; j++)
	  t[i][j]  = 0;
      }

      ZZ_pX T;
      T.rep.SetLength(n);
      for (long j = 0; j < n; j++){
	for (long i = 0; i < nprimes; i++)
	  tmp[i] = t[i][j];
	ZZ_p cT;
	from_modular_rep(cT, tmp, FFTInfo, TmpSpace);
	T.rep[j] = cT;
      }
      T.normalize();
      C2[u][v] = T;
    }
  
  conv(ZC, C2);

  lctx.restore();
  ctx.restore();

  conv(C, ZC);
}
