#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>

#include "magma_output.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* magma output, without assign, using variable var           */
/*------------------------------------------------------------*/
void magma_output(const Mat<zz_pX>& a, const string & var){
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
void magma_assign(const Mat<zz_pX>& a, const string & var, const string & name){
  cout << name << ":=";
  magma_output(a, var);
  cout << ";\n";
}

/*------------------------------------------------------------*/
/* random matrix of a given degree                            */
/*------------------------------------------------------------*/
void random_mat_zz_pX(Mat<zz_pX>& a, long n, long m, long d){
  a.SetDims(n, m);
  for (long i = 0; i < n; i++)
    for (long j = 0; j < m; j++)
      a[i][j] = random_zz_pX(d);
}

/*------------------------------------------------------------*/
/* maximum degree of the entries of a                         */
/*------------------------------------------------------------*/
long deg(const Mat<zz_pX> & a){
  long d = -1;
  for (long i = 0; i < a.NumRows(); i++)
    for (long j = 0; j < a.NumCols(); j++)
      d = max(d, deg(a[i][j]));

  return d;
}


/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* Waksman's algorithm                                        */
/*------------------------------------------------------------*/
void multiply_waksman(Mat<zz_pX> &C, const Mat<zz_pX> &A, const Mat<zz_pX> &B){
  Vec<zz_pX> Arow, Acol, D, E;
  zz_pX val0, val1, val2, crow;

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

  if ( (n&1) == 1){
    for (long l=0; l<m; l++)
      for (long k=0; k<p; k++)
  	C[l][k] += A[l][n-1]*B[n-1][k];
  }
}



/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* naive algorithm                                            */
/*------------------------------------------------------------*/
void multiply_naive(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b){

  long u = a.NumRows();
  long v = a.NumCols();
  long w = b.NumCols();

  c.SetDims(u, w);

  for (long i = 0; i < u; i++)
    for (long j = 0; j < w; j++){
      c[i][j] = 0;
      for (long k = 0; k < v; k++)
	c[i][j] += a[i][k]*b[k][j];
    }
}

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/*------------------------------------------------------------*/
void multiply_evaluate(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b){
  
  long u = a.NumRows();
  long v = a.NumCols();
  long w = b.NumCols();

  c.SetDims(u, w);

  long dA = deg(a);
  long dB = deg(b);
  long dC = dA+dB;
  long sz = dC+1;

  zz_pX_Multipoint * ev;
  zz_pX_Multipoint_Geometric ev_geom;
  zz_pX_Multipoint_FFT ev_FFT;

 if (zz_pInfo->p_info == NULL){ // not an FFT prime
    zz_p aq;
    element_of_order(aq, 2*sz);
    ev_geom = zz_pX_Multipoint_Geometric(aq*aq, sz, 0);
    ev = &ev_geom;
  }
  else{ // FFT prime: always better to use FFT, even if we lose a bit on the number of evaluations
    long k = NextPowerOfTwo(sz);
    ev_FFT = zz_pX_Multipoint_FFT(1L << k);
    ev = &ev_FFT;
  }

  long nb = ev->length();

  Vec<Mat<zz_p>> vala, valb, valc;
  vala.SetLength(nb);
  valb.SetLength(nb);
  valc.SetLength(nb);

  for (long i = 0; i < nb; i++){
    vala[i].SetDims(u, v);
    valb[i].SetDims(v, w);
    valc[i].SetDims(u, w);
  }

  Vec<zz_p> tmp;
  tmp.SetLength(nb);

  // evaluations of A
  for (long i = 0; i < u; i++)
    for (long j = 0; j < v; j++){
      ev->evaluate(tmp, a[i][j]);
      for (long ell = 0; ell < nb; ell++)
	vala[ell][i][j] = tmp[ell];
    }

  // evaluations of B
  for (long i = 0; i < v; i++)
    for (long j = 0; j < w; j++){
      ev->evaluate(tmp, b[i][j]);
      for (long ell = 0; ell < nb; ell++)
	valb[ell][i][j] = tmp[ell];
    }

  // pairwise muls
  for (long ell = 0; ell < nb; ell++)
    mul(valc[ell], vala[ell], valb[ell]);

  // interpolation of C
  for (long i = 0; i < u; i++)
    for (long j = 0; j < w; j++){
      for (long ell = 0; ell < nb; ell++)
	tmp[ell] = valc[ell][i][j];
      ev->interpolate(c[i][j], tmp);
    }
  
}
