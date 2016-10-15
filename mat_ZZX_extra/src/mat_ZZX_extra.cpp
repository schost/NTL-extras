#include <NTL/new.h>
#include <NTL/vec_long.h>
#include <NTL/vec_ulong.h>
#include <NTL/vec_double.h>
#include <gmp.h>
#include "util.h"
#include "vec_vec_vec_zz_pX.h"
#include "mat_lzz_pX_mul.h"
#include "mat_ZZX.h"

// #define CEIL(a,b) (((a)+(b)-1)/(b))
// #define ALLOC(p) (((long *) (p))[0])
// #define SIZE(p) (((long *) (p))[1])
// #define DATA(p) ((mp_limb_t *) (((long *) (p)) + 2))

static inline mp_limb_t * DATA(_ntl_gbigint p) { 
  return ((mp_limb_t *) (((long *) (p)) + 2)); 
}


#define GET_SIZE_NEG(sz, neg, p)  \
do  \
{   \
   long _s;   \
   _s = SIZE(p);   \
   if (_s < 0) {  \
      sz = -_s;  \
      neg = 1;  \
   }  \
   else {  \
      sz = _s;  \
      neg = 0;  \
   }  \
}  \
while (0)

#define STRIP(sz, p)  \
do  \
{  \
   long _i;  \
   _i = sz - 1;  \
   while (_i >= 0 && p[_i] == 0) _i--;  \
   sz = _i + 1;  \
}  \
while (0) 

void ForceNormal(_ntl_gbigint x)
{
   long sx, xneg;
   mp_limb_t *xdata;

   if (!x) return;
   GET_SIZE_NEG(sx, xneg, x);
   xdata = DATA(x);
   STRIP(sx, xdata);
   if (xneg) sx = -sx;
   SIZE(x) = sx;
}


NTL_START_IMPL

NTL_matrix_impl(ZZX,vec_ZZX,vec_vec_ZZX,mat_ZZX)
NTL_io_matrix_impl(ZZX,vec_ZZX,vec_vec_ZZX,mat_ZZX)
NTL_eq_matrix_impl(ZZX,vec_ZZX,vec_vec_ZZX,mat_ZZX)


long deg(const mat_ZZX& a){
  long d = 0;
  long m = a.NumRows();
  long n = a.NumCols();

  for (long i = 0; i < m; i++)
    for (long j = 0; j < n; j++)
      d = max(deg(a[i][j]), d);
  return d;
}

long MaxBits(const mat_ZZX& a){
  long d = 0;
  long m = a.NumRows();
  long n = a.NumCols();

  for (long i = 0; i < m; i++)
    for (long j = 0; j < n; j++)
      d = max(MaxBits(a[i][j]), d);
  return d;
}


/*------------------------------------------------------------*/
/*  direct multiplication                                     */
/*------------------------------------------------------------*/
void mul_direct(mat_ZZX& X, const mat_ZZX& A, const mat_ZZX& B){  
   long n = A.NumRows();  
   long l = A.NumCols();  
   long m = B.NumCols();  
  
   if (l != B.NumRows())  
      Error("matrix mul: dimension mismatch");  
  
   X.SetDims(n, m);  
  
   long i, j, k;  
   ZZX acc, tmp;  
  
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
void mul_waksman(mat_ZZX &C, const mat_ZZX &A, const mat_ZZX &B){
  vec_ZZX Arow, Acol, D, E;
  ZZX val0, val1, val2, crow;

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
  if (n&1 == 1){
    for (long l=0; l<m; l++)
      for (long k=0; k<p; k++)
  	C[l][k] += A[l][n-1]*B[n-1][k];
  }
}

// p1, p2, p3 must hav been initialized with the proper size = (deg(F)+1)
void to_CRT(vec_vec_vec_zz_pX& p1, vec_vec_vec_zz_pX& p2, vec_vec_vec_zz_pX& p3, const mat_ZZX& F,
	    const long q1, const long q2, const long q3){
  long d = deg(F);
  long h = CEIL(MaxBits(F), NTL_BITS_PER_LONG);

  long m = F.NumRows();
  long n = F.NumCols();

  p1.SetLength(h);
  p2.SetLength(h);
  p3.SetLength(h);

  for (long i = 0; i < h; i++){
    p1[i].SetLength(m);
    p2[i].SetLength(m);
    p3[i].SetLength(m);
    for (long j = 0; j < m; j++){
      p1[i][j].SetLength(n);
      p2[i][j].SetLength(n);
      p3[i][j].SetLength(n);
      for (long k = 0; k < n; k++){
	p1[i][j][k].rep.SetLength(d+1);
	p2[i][j][k].rep.SetLength(d+1);
	p3[i][j][k].rep.SetLength(d+1);
      }
    }
  }

  for (long a = 0; a < m; a++)
    for (long b = 0; b < n; b++)
      for (long i = 0; i <= d; i++){
  	ZZ tmp = coeff(F[a][b], i);
  	if (tmp == 0)
  	  for (long j = 0; j < h; j++){
  	    p1[j][a][b].rep[i]._zz_p__rep = 0;
  	    p2[j][a][b].rep[i]._zz_p__rep = 0;
  	    p3[j][a][b].rep[i]._zz_p__rep = 0;
  	  }
  	else{
  	  mp_limb_t *dat = DATA(tmp.rep);
  	  long sz = SIZE(tmp.rep);
	  // todo: check sign
  	  for (long j = 0; j < sz; j++){
  	    p1[j][a][b].rep[i]._zz_p__rep = dat[j] % q1;
  	    p2[j][a][b].rep[i]._zz_p__rep = dat[j] % q2;
  	    p3[j][a][b].rep[i]._zz_p__rep = dat[j] % q3;
  	  }
  	  for (long j = sz; j < h; j++){
  	    p1[j][a][b].rep[i]._zz_p__rep = 0;
  	    p2[j][a][b].rep[i]._zz_p__rep = 0;
  	    p3[j][a][b].rep[i]._zz_p__rep = 0;
  	  }
  	}
      }
}


void from_CRT(mat_ZZX& C, const vec_vec_vec_zz_pX& p1, const vec_vec_vec_zz_pX& p2, const vec_vec_vec_zz_pX& p3, 
	      const long d, const long q1, const long q2, const long q3){

  long h = p1.length();
  long m = p1[0].length();
  long n = p1[0][0].length();

  ZZ zp1 = to_ZZ(q1), zp2 = to_ZZ(q2), zp3 = to_ZZ(q3);
  ZZ mm = zp1*zp2*zp3;
  ZZ zP1 = zp2*zp3*InvMod((zp2*zp3) % zp1, zp1);
  ZZ zP2 = zp1*zp3*InvMod((zp1*zp3) % zp2, zp2);
  ZZ zP3 = zp1*zp2*InvMod((zp1*zp2) % zp3, zp3);

  mp_limb_t q[4];
  mp_limb_t res[4];
  mp_limb_t extra;

  mp_limb_t *d1 = DATA(zP1.rep);
  mp_limb_t *d2 = DATA(zP2.rep);
  mp_limb_t *d3 = DATA(zP3.rep);
  mp_limb_t *dm = DATA(mm.rep);

  C.SetDims(m, n);
  for (long a = 0; a < m; a++){
    for (long b = 0; b < n; b++){
      C[a][b].rep.SetLength(d+1);
      for (long i = 0; i <= d; i++){
	ZZ tmp(INIT_SIZE, h+3);
	mp_limb_t* dat = DATA(tmp.rep);
	SIZE(tmp.rep) = h+3;
	for (long j = 0; j < h+3; j++)
	  dat[j] = 0;
	for (long j = 0; j < h; j++){
	  extra = mpn_mul_1(res, d1, 3, p1[j][a][b].rep[i]._zz_p__rep);
	  res[3] = extra;
	  extra = mpn_addmul_1(res, d2, 3, p2[j][a][b].rep[i]._zz_p__rep);
	  res[3] += extra;
	  extra = mpn_addmul_1(res, d3, 3, p3[j][a][b].rep[i]._zz_p__rep);
	  res[3] += extra;
	  mpn_tdiv_qr (q, res, 0L, res, 4, dm, 3);
	  mpn_add_n(dat+j, dat+j, res, 3); // do not need the carry?
	}
	ForceNormal(tmp.rep);
	C[a][b].rep[i] = ZZ(tmp);
      }
      C[a][b].normalize();
    }
  }
}


void mul_CRT(mat_ZZX &C, const mat_ZZX &A, const mat_ZZX &B){
  long m = A.NumRows();
  long n = A.NumCols();
  long p = B.NumCols();

  vec_vec_vec_zz_pX a1, a2, a3;
  vec_vec_vec_zz_pX b1, b2, b3;
  vec_vec_vec_zz_pX c1, c2, c3;


  // precompute all this stuff
  zz_p::FFTInit(0);
  long p1 = zz_p::modulus();
  zz_p::FFTInit(1);
  long p2 = zz_p::modulus();
  zz_p::FFTInit(2);
  long p3 = zz_p::modulus();

  INITTIME;

  GETTIME;
  to_CRT(a1, a2, a3, A, p1, p2, p3);
  to_CRT(b1, b2, b3, B, p1, p2, p3);
  PRINTIDX(" CRT    ");
  PRINTTIME;


  GETTIME;
  zz_p::FFTInit(0);
  long h1 = b1.length();
  c1.SetLength(h1);

  for (long a = 0; a < m; a++)
    for (long b = 0; b < n; b++)
      a1[0][a][b].normalize();
    
  for (long i = 0; i < h1; i++){
    for (long a = 0; a < n; a++)
      for (long b = 0; b < p; b++)
	b1[i][a][b].normalize();

    multiply_FFT(c1[i], a1[0], b1[i]);
  }    

  zz_p::FFTInit(1);
  long h2 = b2.length();
  c2.SetLength(h2);

  for (long a = 0; a < m; a++)
    for (long b = 0; b < n; b++)
      a2[0][a][b].normalize();
    
  for (long i = 0; i < h2; i++){
    for (long a = 0; a < n; a++)
      for (long b = 0; b < p; b++)
	b2[i][a][b].normalize();

    multiply_FFT(c2[i], a2[0], b2[i]);
  }    

  zz_p::FFTInit(2);
  long h3 = b3.length();
  c3.SetLength(h3);

  for (long a = 0; a < m; a++)
    for (long b = 0; b < n; b++)
      a3[0][a][b].normalize();
    
  for (long i = 0; i < h3; i++){
    for (long a = 0; a < n; a++)
      for (long b = 0; b < p; b++)
	b3[i][a][b].normalize();

    multiply_FFT(c3[i], a3[0], b3[i]);
  }    
  PRINTTIME;

  GETTIME;
  from_CRT(C, c1, c2, c3, deg(A)+deg(B), p1, p2, p3);
  PRINTTIME;
  PRINTNL;
}


NTL_END_IMPL
