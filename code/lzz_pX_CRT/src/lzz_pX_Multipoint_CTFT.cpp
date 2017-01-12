#include "lzz_pX_CRT.h"

NTL_CLIENT

/*-----------------------------------------------------------*/
/*-----------------------------------------------------------*/
/* init                                                      */
/*-----------------------------------------------------------*/
/*-----------------------------------------------------------*/

/*-----------------------------------------*/
/* returns the exponents ki such that      */
/*    n = sum_i 2^k_i                      */
/*-----------------------------------------*/
void CTFT_exponents(Vec<long>& vec, const long nn){

  long i = weight(nn);
  vec.SetLength(i+1);  

  if (nn == 0){
    vec[0] = -1;
    return;
  }
  
  long j;
  for (j = 0; j < i; j++)
    vec[j] = 0;
  vec[i] = -1;

  long n = nn;
  long a = 1;
  i = 0;
  j = 0;
  while (a <= n/2){
    a <<= 1;
    j++;
  }
  vec[i++] = j;

  while (a != n){
    long k = j-1;
    long b = a >> 1;
    while (! (n & b)){
      k--;
      b >>= 1;
    }
    n -= a;
    a = b;
    vec[i++] = k;
    j = k;
  }
}

/*-----------------------------------------------------------*/
/* sets all entries in the vector pts                        */ 
/*-----------------------------------------------------------*/
void zz_pX_Multipoint_CTFT::init_points(){
  Vec<long> degrees;
  CTFT_exponents(degrees, n);
  pts.SetLength(n);

  long j = 0;
  long acc = 0;

  zz_pInfoT *info = zz_pInfo;
  FFTPrimeInfo* p_info = info->p_info;
  const long * root = p_info->RootTable[0].elts();

  while (degrees[j] != -1){
    zz_p r0 = to_zz_p(root[degrees[j]]);
    zz_p r1 = to_zz_p(root[degrees[j]+1]);
    zz_p r2 = to_zz_p(1);

    for (long i = 0; i < (1L << degrees[j]); i++){
      pts[acc+i] = r1*r2;
      r2 *= r0;
    }

    acc += (1L << degrees[j]);
    j++;
  }
}

/*--------------------------------------------------------------*/
/* initializes the k-th row of pre-multipliers for negacyclic   */
/* convolution (== in size 2^k)                                 */
/*--------------------------------------------------------------*/
void zz_pX_Multipoint_CTFT::init_multipliers(const long k){

   if (k <= MaxK) 
     return;

   z.SetLength(k+1);
   z_precomp.SetLength(k+1);
   invz.SetLength(k+1);
   invz_precomp.SetLength(k+1);

   // TODO: check that non NULL
   zz_pInfoT *info;
   FFTPrimeInfo* p_info;
   const long * root;
   info = zz_pInfo;
   p_info = info->p_info;
   root = p_info->RootTable[0].elts();

   long p = zz_p::modulus();
   mulmod_t pinv = zz_p::ModulusInverse();

   // root[i] is a root of order 2^i
   // so root[0] = 1, root[1] = -1
   //
   // root[i+1]^j, j = 0..2^i-1
   for (long i = MaxK+1; i <= k; i++){
     z[i].SetLength(1 << i);
     z_precomp[i].SetLength(1 << i);
     invz[i].SetLength(1 << i);
     invz_precomp[i].SetLength(1 << i);

     zz_p tmpz = to_zz_p(root[i+1]);
     zz_p tmp_invz = 1/tmpz;

     z[i][0] = 1;
     z_precomp[i][0] = PrepMulModPrecon(z[i][0], p, pinv);
     invz[i][0] = p_info->TwoInvTable[i];
     invz_precomp[i][0] = PrepMulModPrecon(invz[i][0], p, pinv);

     for (long j = 1; j < (1L << i); j++){
       z[i][j] = (to_zz_p(z[i][j-1]) * tmpz).LoopHole();
       z_precomp[i][j] = PrepMulModPrecon(z[i][j], p, pinv);
       invz[i][j] = (to_zz_p(invz[i][j-1]) * tmp_invz).LoopHole();
       invz_precomp[i][j] = PrepMulModPrecon(invz[i][j], p, pinv);
     }
   }

   MaxK = k;
}


/*-----------------------------------------------------------*/
/*-----------------------------------------------------------*/
/* forward algorithms                                        */
/*-----------------------------------------------------------*/
/*-----------------------------------------------------------*/

/*-----------------------------------------------------------*/
/* applies the multi-reduction map to x (mod p), in place    */
/* input is in [0,p), output is in [0,2p)                    */
/* x must have size at least 2N, and upper half must be zero */
/*-----------------------------------------------------------*/
void zz_pX_Multipoint_CTFT::reduce(long* x) const {

  const long p = zz_p::modulus();
  const long N = n;
  if (N <= 1)
    return;

  long nn = N;
  long k = 0;
  long a = 1L;
  while (a <= nn/2){
    k++;
    a <<= 1;
  }

  // if N is a power of 2, nothing to do
  if (a == nn)
    return;

  // else, loop through the reduction process
  while (a != nn){

    long b = a >> 1;
    long ell = k-1;
    while (! (nn & b)){
      b >>= 1;
      ell--;
    }

    long t = 0;
    const long b2 = b << 1;

    // b = 1: unroll the loop
    if (b == 1){
      long u, v;
      u = x[0];
      v = x[a+0];
      x[0] = u - v + p;
      x[a+0] = AddMod(u, v, p);
      u = x[1];
      v = x[a+1];
      x[1] = u - v + p;
      x[a+1] = AddMod(u, v, p);

      long t = 2;
      long lambda = 1L << (k-ell-1);
      for (long i = 1; i < lambda; i++){
	long u, v;
	u = x[t];
	v = x[a+t];
	x[t] = u - v + p;
	x[a] = AddMod(AddMod(u, v, p), x[a], p);
	t++;
	u = x[t];
	v = x[a+t];
	x[t] = u - v + p;
	x[a+1] = AddMod(AddMod(u, v, p), x[a+1], p);
	t++;
      }
    }
    else{
      // i = 0 to 2^(k-ell-1) = a/(2b), j = 0 to 2b
      // i = 0 : special case
      for (; t < b2; t++){
	long u, v;
	u = x[t];
	v = x[a+t];
	x[t] = u - v + p;
	x[a+t] = AddMod(u, v, p);
      }
      
      long lambda = 1L << (k-ell-1);
      for (long i = 1; i < lambda; i++){
	for (long j = 0; j < b2; j++){
	  long u, v;
	  u = x[t];
	  v = x[a+t];
	  x[t] = u - v + p;
	  x[a+j] = AddMod(AddMod(u, v, p), x[a+j], p);
	  t++;
	}
      }
    }

    x += a;
    nn -= a;
    a = b;
    k = ell;
  }

  for (long t = 0; t < a; t++){
    long u, v;
    u = x[t];
    v = x[a+t];
    x[t] = u - v + p;
  }
}

/*---------------------------------------------------------*/
/* reduces s mod (X^sz-1), assuming s has length N         */
/* assumes sz divides N                                    */
/* all calculations are mod p                              */
/*---------------------------------------------------------*/
static inline 
void fold_minus(long *x, const long *s, const long sz, const long N, const long p){
  long i;

  for (i = 0; i < sz; i++)
    x[i] = s[i];

  if (sz == 2){
    while (i < N){
      x[0] = AddMod(x[0], s[i], p);
      i++;
      x[1] = AddMod(x[1], s[i], p);
      i++;
    }
  }
  else{
    while (i < N)
      for (long j = 0; j < sz; j++, i++)
	x[j] = AddMod(x[j], s[i], p);
  }
}

/*-----------------------------------------------------------*/
/* applies the CRT map to x (mod p), in place                */
/* x has length n                                            */
/* result is in y                                            */
/* input is in [0,p), output is in [0,p)                     */
/* tmp is a temporary workspace, of length at least 2n       */
/*-----------------------------------------------------------*/
void zz_pX_Multipoint_CTFT::CRT(long* x, long *tmp) const{

  const long p = zz_p::modulus();
  const mulmod_t pinv = zz_p::ModulusInverse();
  const long N = n;

  // p must be odd!
  // half = 1/2 mod p
  const long half = NegateMod(p >> 1, p);
  const mulmod_precon_t half_pinv = PrepMulModPrecon(half, p, pinv);

  if (N <= 1)
    return;

  long nn = N;
  long a, b, b2, n2;
  a = 1L;
  n2 = (nn >> 1);
  while (a <= n2)
    a <<= 1;

  while (nn != a){
    nn = nn-a;
    b = 1L;
    n2 = (nn >> 1);
    while (b <= n2)
      b <<= 1;
    b2 = b << 1;
    
    fold_minus(tmp, x, b2, a, p);

    for(long i = 0; i < b; i++)
      x[a+i] = AddMod(AddMod(tmp[b+i], tmp[b+i], p), x[a+i], p);

    tmp += b2;
    x += a;
    a = b;
  }

  while (nn != N){
    b = a;
    a <<= 1;
    b2 = a;
    while (!(a & N))
      a <<= 1;
    tmp -= b2;
    x -= a;
    
    for (long i = 0; i < b; i++){
      long u = AddMod(tmp[i], tmp[b+i], p);
      u = SubMod(x[a+i], u, p);
      x[a+i] = MulModPrecon(u, half, p, half_pinv);
      x[i] = AddMod(x[i], x[a+i], p);
    }
      
    for (long i = b; i < nn; i++){
      x[a+i] = MulModPrecon(x[a+i], half, p, half_pinv);
      x[i] = AddMod(x[i], x[a+i], p);
    }

    nn = nn+a;

  }
}

/*-----------------------------------------------------------*/
/* evaluates a chunk of length 2^k                           */
/* applies the premultiplier before FFT                      */
/* input length = 2^k                                        */
/* input is in [0,2p), output is in [0,p)                    */
/*-----------------------------------------------------------*/
void zz_pX_Multipoint_CTFT::evaluate_chunk(zz_p * val, const long * coeffs, const long k, long * wk) const{
  
  long K = 1L << k;

  if (k == 0){
    val[0] = coeffs[0];
    return;
  }
  
  const long p = zz_p::modulus();
  zz_pInfoT *info = zz_pInfo;
  FFTPrimeInfo *p_info = info->p_info;

  const long * powers = z[k].elts();
  const mulmod_precon_t * powers_precomp = z_precomp[k].elts();

  for (long i = 0; i < K; i++){
    long v = coeffs[i] > p ? coeffs[i] - p : coeffs[i];
    wk[i] = MulModPrecon(v, powers[i], p, powers_precomp[i]);
  }

  FFTFwd(wk, wk, k, *p_info);
  
  for (long i = 0; i < K; i++)
    val[i] = to_zz_p(wk[i]);
}


/*-----------------------------------------------------------*/
/* interpolates a chunk of length 2^k                        */
/* input length = 2^k                                        */
/* input is word-size, output is in [0,p)                    */
/*-----------------------------------------------------------*/
void zz_pX_Multipoint_CTFT::interpolate_chunk(long * val, const zz_p * coeffs, const long k, long * wk) const{
  
  long K = 1L << k;

  if (k == 0){
    val[0] = coeffs[0]._zz_p__rep;
    return;
  }
  
  const long p = zz_p::modulus();
  zz_pInfoT *info = zz_pInfo;
  FFTPrimeInfo *p_info = info->p_info;

  const long * powers = invz[k].elts();
  const mulmod_precon_t * powers_precomp = invz_precomp[k].elts();

  for (long i = 0; i < K; i++)
    wk[i] = coeffs[i]._zz_p__rep;
  
  FFTRev(wk, wk, k, *p_info);

  for (long i = 0; i < K; i++)
    val[i] = MulModPrecon(wk[i], powers[i], p, powers_precomp[i]);
}


/*---------------------------------------------------------*/
/* evaluates A at n points                                 */
/*---------------------------------------------------------*/
void zz_pX_Multipoint_CTFT::evaluate(Vec<zz_p>& val, const long * A, const long n0) const {

  val.SetLength(n);

  if (n == 0)
    return;

  if (n == 1){
    val[0] = A[0];
    return;
  }

  zz_p * val_elts = val.elts();

  Vec<long> wk_vec;
  wk_vec.SetLength(3*n);
  long *wk = wk_vec.elts();
  long *wk2 = wk + 2*n;

  for (long i = 0; i < n0; i++)
    wk[i] = A[i];
  for (long i = n0; i < 3*n; i++)
    wk[i] = 0;
  reduce(wk);

  long nn = n;
  long k = 0;
  long aa = 1L;

  while (aa <= nn/2){
    k++;
    aa <<= 1;
  }

  const long p = zz_p::modulus();
  zz_pInfoT *info = zz_pInfo;
  FFTPrimeInfo *p_info = info->p_info;

  do {
    if (k == 0)
      val_elts[0] = wk[0];
    else {
      const long * powers = z[k].elts();
      const mulmod_precon_t * powers_precomp = z_precomp[k].elts();
      for (long i = 0; i < aa; i++){
	wk2[i] = sp_CorrectExcess(MulModPrecon(wk[i], powers[i], p, powers_precomp[i]), p);
	// long v = wk[i] > p ? wk[i] - p : wk[i];
	// wk2[i] = MulModPrecon(v, powers[i], p, powers_precomp[i]);
      }
      FFTFwd(wk2, wk2, k, *p_info);
      for (long i = 0; i < aa; i++)
	val_elts[i]._zz_p__rep = wk2[i];
    }
    
    val_elts += aa;
    wk += aa;

    nn = nn-aa;

    aa = 1L;
    k = 0;
    while (aa <= nn/2){
      k++;
      aa <<= 1;
    }
  }
  while (nn != 0);
}

/*---------------------------------------------------------*/
/* evaluates A at n points                                 */
/*---------------------------------------------------------*/
void zz_pX_Multipoint_CTFT::evaluate(Vec<zz_p>& val, const zz_pX& A) const {

  const zz_p *a = A.rep.elts();
  Vec<long> vA;
  long n0 = A.rep.length();
  vA.SetLength(n0);
  long * vAe = vA.elts();
  for (long i = 0; i < n0; i++)
    vAe[i] = a[i]._zz_p__rep;

  evaluate(val, vAe, n0);
}

/*---------------------------------------------------------*/
/* interpolates A at n points                              */
/*---------------------------------------------------------*/
void zz_pX_Multipoint_CTFT::interpolate(zz_pX& f, const Vec<zz_p>& val) const{
  if (n == 0){
    f = 0;
    return;
  }

  if (n == 1){
    SetCoeff(f, 0, val[0]);
    return;
  }

  const zz_p * val_elts = val.elts();

  Vec<long> wk_vec;
  wk_vec.SetLength(3*n);
  long *wk = wk_vec.elts();
  long *wk2 = wk + n;
  for (long i = 0; i < 3*n; i++)
    wk[i] = 0;
  
  long nn = n;
  long k = 0;
  long aa = 1L;

  const long p = zz_p::modulus();
  zz_pInfoT *info = zz_pInfo;
  FFTPrimeInfo *p_info = info->p_info;

  while (aa <= nn/2){
    k++;
    aa <<= 1;
  }

  do {
    if (k == 0){
      wk[0] = val_elts[0]._zz_p__rep;
    }
    else {
      const long * powers = invz[k].elts();
      const mulmod_precon_t * powers_precomp = invz_precomp[k].elts();

      for (long i = 0; i < aa; i++)
	wk2[i] = val_elts[i]._zz_p__rep;

      FFTRev(wk2, wk2, k, *p_info);

      for (long i = 0; i < aa; i++)
	wk[i] = MulModPrecon(wk2[i], powers[i], p, powers_precomp[i]);
    }
    val_elts += aa;
    wk += aa;
    nn = nn-aa;

    aa = 1L;
    k = 0;
    while (aa <= nn/2){
      k++;
      aa <<= 1;
    }
  }
  while (nn != 0);

  CRT(wk_vec.elts(), wk2);

  f.rep.SetLength(n);
  for (long i = 0; i < n; i++)
    f.rep[i] = to_zz_p(wk_vec[i]);
  f.normalize();
}

/*---------------------------------------------------------*/
/*---------------------------------------------------------*/
/* TFT multiplication                                      */
/* same loop for all x, y and z; saves a bit               */
/*---------------------------------------------------------*/
/*---------------------------------------------------------*/

/*-----------------------------------------------------------*/
/* applies the multi-reduction map to x,y (mod p), in place  */
/* input is in [0,p), output is in [0,2p)                    */
/* x must have size at least 2N, and upper half must be zero */
/*-----------------------------------------------------------*/
void zz_pX_Multipoint_CTFT::bi_reduce(long* x, long* y) const {

  const long p = zz_p::modulus();
  const long N = n;
  if (N <= 1)
    return;

  long nn = N;
  long k = 0;
  long a = 1L;
  while (a <= nn/2){
    k++;
    a <<= 1;
  }

  // if N is a power of 2, nothing to do
  if (a == nn)
    return;

  // else, loop through the reduction process
  while (a != nn){

    long b = a >> 1;
    long ell = k-1;
    while (! (nn & b)){
      b >>= 1;
      ell--;
    }

    long t = 0;
    const long b2 = b << 1;

    // b = 1: unroll the loop
    if (b == 1){
      long u, v;
      u = x[0];
      v = x[a+0];
      x[0] = u - v + p;
      x[a+0] = AddMod(u, v, p);
      u = x[1];
      v = x[a+1];
      x[1] = u - v + p;
      x[a+1] = AddMod(u, v, p);
      u = y[0];
      v = y[a+0];
      y[0] = u - v + p;
      y[a+0] = AddMod(u, v, p);
      u = y[1];
      v = y[a+1];
      y[1] = u - v + p;
      y[a+1] = AddMod(u, v, p);

      long t = 2;
      long lambda = 1L << (k-ell-1);
      for (long i = 1; i < lambda; i++){
	long u, v;
	u = x[t];
	v = x[a+t];
	x[t] = u - v + p;
	x[a] = AddMod(AddMod(u, v, p), x[a], p);
	u = y[t];
	v = y[a+t];
	y[t] = u - v + p;
	y[a] = AddMod(AddMod(u, v, p), y[a], p);
	t++;
	u = x[t];
	v = x[a+t];
	x[t] = u - v + p;
	x[a+1] = AddMod(AddMod(u, v, p), x[a+1], p);
	u = y[t];
	v = y[a+t];
	y[t] = u - v + p;
	y[a+1] = AddMod(AddMod(u, v, p), y[a+1], p);
	t++;
      }
    }
    else{
      // i = 0 to 2^(k-ell-1) = a/(2b), j = 0 to 2b
      // i = 0 : special case
      for (; t < b2; t++){
	long u, v;
	u = x[t];
	v = x[a+t];
	x[t] = u - v + p;
	x[a+t] = AddMod(u, v, p);
	u = y[t];
	v = y[a+t];
	y[t] = u - v + p;
	y[a+t] = AddMod(u, v, p);
      }
      
      long lambda = 1L << (k-ell-1);
      for (long i = 1; i < lambda; i++){
	for (long j = 0; j < b2; j++){
	  long u, v;
	  u = x[t];
	  v = x[a+t];
	  x[t] = u - v + p;
	  x[a+j] = AddMod(AddMod(u, v, p), x[a+j], p);
	  u = y[t];
	  v = y[a+t];
	  y[t] = u - v + p;
	  y[a+j] = AddMod(AddMod(u, v, p), y[a+j], p);
	  t++;
	}
      }
    }

    x += a;
    y += a;
    nn -= a;
    a = b;
    k = ell;
  }

  for (long t = 0; t < a; t++){
    long u, v;
    u = x[t];
    v = x[a+t];
    x[t] = u - v + p;
    u = y[t];
    v = y[a+t];
    y[t] = u - v + p;
  }
}

/*---------------------------------------------------------*/
/* negacyclic convolution in length 2^k                    */
/* wk and wk2 are workspaces of length 2^k                 */
/*---------------------------------------------------------*/
void zz_pX_Multipoint_CTFT::negacyclic_convolution(long *c, const long *a, const long *b, const long k, long *wk, long *wk2) const{
  
  const long p = zz_p::modulus();

  switch(k){
  case 0:
    c[0] = MulMod(a[0], b[0], p);
    return;
  case 1:
    c[0] = SubMod(MulMod(a[0], b[0], p), MulMod(a[1], b[1], p), p);
    c[1] = AddMod(MulMod(a[0], b[1], p), MulMod(a[1], b[0], p), p);
    return;
  default:
    zz_pInfoT *info = zz_pInfo;
    FFTPrimeInfo *p_info = info->p_info;
    
    const long * powers = z[k].elts();
    const mulmod_precon_t * powers_precomp = z_precomp[k].elts();
    const long * inv_powers = invz[k].elts();
    const mulmod_precon_t * inv_powers_precomp = invz_precomp[k].elts();
    
    long K = 1L << k;

    const mulmod_t pinv = zz_p::ModulusInverse();

    for (long i = 0; i < K; i++){
      long v;
      v = a[i] > p ? a[i] - p : a[i];
      wk[i] = MulModPrecon(v, powers[i], p, powers_precomp[i]);
      v = b[i] > p ? b[i] - p : b[i];
      wk2[i] = MulModPrecon(v, powers[i], p, powers_precomp[i]);
    }

    FFTFwd(wk, wk, k, *p_info);
    FFTFwd(wk2, wk2, k, *p_info);
    
    for (long i = 0; i < K; i++)
      c[i] = MulMod(wk[i], wk2[i], p, pinv);
    
    FFTRev(c, c, k, *p_info);

    for (long i = 0; i < K; i++)
      c[i] = MulModPrecon(c[i], inv_powers[i], p, inv_powers_precomp[i]);
  }
}

/*---------------------------------------------------------*/
/* multiplies a and b                                      */
/*---------------------------------------------------------*/
void zz_pX_Multipoint_CTFT::mul(zz_pX& C, const zz_pX& A, const zz_pX& B) const {
  if (A == 0 || B == 0){
    C = 0;
    return;
  }

  long nA = deg(A) + 1;
  long nB = deg(B) + 1;
  const long N = n;
  long p = zz_p::modulus();

  Vec<long> wk_vec;
  wk_vec.SetLength(5*N);
  long *wk =wk_vec.elts();
  long *wkN = wk + N;
  long *wkN2 = wk + 2*N;
  long *wkN3 = wk + 3*N;
  long *wkN4 = wk + 4*N;
  long *wk_bak = wk;

  for (long i = 0; i < 5*N; i++)
    wk[i] = 0;

  const zz_p *a = A.rep.elts();
  for (long i = 0; i < nA; i++)
    wk[i] = a[i]._zz_p__rep;

  const zz_p *b = B.rep.elts();
  for (long i = 0; i < nB; i++)
    wkN2[i] = b[i]._zz_p__rep;
  bi_reduce(wk, wkN2);
  // reduce(wk);
  // reduce(wkN2);


  if (N == 1)
    wkN4[0] = MulMod(wk[0], wkN2[0], p);
  else {
    long nloc = N;
    long k = 0;
    long a = 1L;
    
    while (a <= nloc/2){
      k++;
      a <<= 1;
    }
    
    do {
      negacyclic_convolution(wkN4, wk, wkN2, k, wkN, wkN3);
      
      wkN4 += a;
      wk += a;
      wkN2 += a;
      
      nloc = nloc-a;
      
      a = 1L;
      k = 0;
      while (a <= nloc/2){
	k++;
	a <<= 1;
      }
    }
    while (nloc != 0);
  }

  wk = wk_bak;
  wkN4 = wk + 4*N;
  
  CRT(wkN4, wk);

  clear(C);
  C.rep.SetLength(N);
  zz_p *c = C.rep.elts();
  for (long i = 0; i < N; i++)
    c[i]._zz_p__rep = wkN4[i];
  C.normalize();

}

/*-----------------------------------------------------------*/
/*-----------------------------------------------------------*/
/* transposed algorithms                                     */
/*-----------------------------------------------------------*/
/*-----------------------------------------------------------*/

/*-----------------------------------------------------------*/
/* transpose multi-reduction (= extension of sequences)      */
/* in place                                                  */
/* input has length N, but needs space 2N                    */
/*-----------------------------------------------------------*/
void zz_pX_Multipoint_CTFT::reduce_t(long* x) const {

  const long p = zz_p::modulus();
  const long N = n;

  if (N <= 1)
    return;

  long a = 1;
  long k = 0;

  while (! (a & N)){
    k++;
    a <<= 1;
  }

  if (a == N)
    return;

  for (long i = 0; i < a; i++)
    x[N+i] = NegateMod(x[N-a+i], p);

  long n = a;
  while (n != N){
    long b = a;
    long ell = k;
    a <<= 1;
    k++;
    while (!(a & N)){
      a <<= 1;
      k++;
    }

    const long b2 = b << 1;
    long lambda = 1L << (k-ell-1);

    n = n+a;

    long t = b2;
    // we have to deal with i=0 separately
    // (since otherwise we would erase the source)
    //
    // we also write special cases for b=1 (and should do b=2)
    if (b == 1){
      for (long i = 1; i < lambda; i++){
  	long u, v;
  	u = x[N-n+t];
  	v = x[N-n+a];
  	x[N-n+t] = AddMod(v, u, p);
  	x[N-n+a+t] = SubMod(v, u, p);
  	t++;
  	u = x[N-n+t];
  	v = x[N-n+a+1];
  	x[N-n+t] = AddMod(v, u, p);
  	x[N-n+a+t] = SubMod(v, u, p);
  	t++;
      }
      
      long u, v;
      u = x[N-n];
      v = x[N-n+a];
      x[N-n] = AddMod(v, u, p);
      x[N-n+a] = SubMod(v, u, p);
      u = x[N-n+1];
      v = x[N-n+a+1];
      x[N-n+1] = AddMod(v, u, p);
      x[N-n+a+1] = SubMod(v, u, p);

    }
    else{
      for (long i = 1; i < lambda; i++)
  	for (long j = 0; j < b2; j++){
  	  long u = x[N-n+t];
  	  long v = x[N-n+a+j];
  	  x[N-n+t] = AddMod(v, u, p);
  	  x[N-n+a+t] = SubMod(v, u, p);
  	  t++;
  	}
      
      for (long j = 0; j < b2; j++){
  	long u = x[N-n+j];
  	long v = x[N-n+a+j];
  	x[N-n+j] = AddMod(v, u, p);
  	x[N-n+a+j] = SubMod(v, u, p);
      }
    }
  }

}

/*-----------------------------------------------------------*/
/* transpose CRT (= decomposition of sequences)              */
/* in place                                                  */
/* input has length N                                        */
/*-----------------------------------------------------------*/
void zz_pX_Multipoint_CTFT::CRT_t(long* x, long *tmp) const {
  
  const long N = n;
  const long p = zz_p::modulus();
  const mulmod_t pinv = zz_p::ModulusInverse();

  // p must be odd!
  // half = 1/2 mod p
  const long half = NegateMod(p >> 1, p);
  const mulmod_precon_t half_pinv = PrepMulModPrecon(half, p, pinv);

  if (N <= 1)
    return;

  long n = N;
  long a, b, b2, n2;
  a = 1L;
  n2 = (n >> 1);
  while (a <= n2)
    a <<= 1;

  while (n != a){
    n = n-a;
    b = 1L;
    n2 = (n >> 1);
    while (b <= n2)
      b <<= 1;
    b2 = b << 1;
    
    for (long i = 0; i < b; i++){
      long u = MulModPrecon(AddMod(x[i], x[a+i], p), half, p, half_pinv);
      x[a+i] = u;
      tmp[i] = u;
    }

    for (long i = b; i < n; i++)
      x[a+i] = MulModPrecon(AddMod(x[i], x[a+i], p), half, p, half_pinv);

    tmp += b2;
    x += a;
    a = b;
  }

  while (n != N){
    b = a;
    a <<= 1;
    b2 = a;
    while (!(a & N))
      a <<= 1;
    tmp -= b2;
    x -= a;
    
    for (long i = 0; i < b; i++)
      tmp[i+b] = SubMod(tmp[i], AddMod(x[a+i], x[a+i], p), p);

    long t = 0;
    while (t < a)
      for (long i = 0; i < b2; i++, t++)
    	x[t] = SubMod(x[t], tmp[i], p);

    n = n+a;
  }
}

/*-----------------------------------------------------------*/
/* transpose-evaluates a chunk of length 2^k                 */
/* applies the premultiplier after FFT                       */
/* input length = 2^k                                        */
/*-----------------------------------------------------------*/
void zz_pX_Multipoint_CTFT::evaluate_chunk_t(zz_p * val, const long * coeffs, const long k, long * wk) const{
  
  long K = 1L << k;

  if (k == 0){
    val[0] = coeffs[0];
    return;
  }
  
  const long p = zz_p::modulus();
  zz_pInfoT *info = zz_pInfo;
  FFTPrimeInfo *p_info = info->p_info;

  const long * powers = z[k].elts();
  const mulmod_precon_t * powers_precomp = z_precomp[k].elts();

  for (long i = 0; i < K; i++)
    wk[i] = coeffs[i];

  FFTFwd(wk, wk, k, *p_info);

  for (long i = 0; i < K; i++)
    val[i] = to_zz_p(MulModPrecon(wk[i], powers[i], p, powers_precomp[i]));
}

/*-----------------------------------------------------------*/
/* transpose-interpolates a chunk of length 2^k              */
/* input length = 2^k                                        */
/*-----------------------------------------------------------*/
void zz_pX_Multipoint_CTFT::interpolate_chunk_t(long * val, const zz_p * coeffs, const long k, long * wk) const{
   
  long K = 1L << k;

  if (k == 0){
    val[0] = coeffs[0]._zz_p__rep;
    return;
  }
  
  const long p = zz_p::modulus();
  zz_pInfoT *info = zz_pInfo;
  FFTPrimeInfo *p_info = info->p_info;

  const long * powers = invz[k].elts();
  const mulmod_precon_t * powers_precomp = invz_precomp[k].elts();

  for (long i = 0; i < K; i++)
    wk[i] = MulModPrecon(coeffs[i]._zz_p__rep, powers[i], p, powers_precomp[i]);
  
  FFTRev(val, wk, k, *p_info);
}

/*---------------------------------------------------------*/
/* transpose-evaluates A at n points                       */
/*---------------------------------------------------------*/
void zz_pX_Multipoint_CTFT::evaluate_t(Vec<zz_p>& val, const Vec<zz_p>& A) const {

  const zz_p *a = A.elts();

  val.SetLength(n);
  zz_p * val_elts = val.elts();

  if (n == 0)
    return;

  if (n == 1){
    val_elts[0] = a[0];
    return;
  }

  Vec<long> wk_vec;
  wk_vec.SetLength(3*n);
  long *wk = wk_vec.elts();
  long *wk2 = wk + 2*n;

  for (long i = 0; i < n; i++)
    wk[i] = a[i]._zz_p__rep;
  for (long i = n; i < 3*n; i++)
    wk[i] = 0;

  long nn = n;
  long k = 0;
  long aa = 1L;

  while (aa <= nn/2){
    k++;
    aa <<= 1;
  }

  do {
    evaluate_chunk_t(val_elts, wk, k, wk2);
    val_elts += aa;
    wk += aa;

    nn = nn-aa;

    aa = 1L;
    k = 0;
    while (aa <= nn/2){
      k++;
      aa <<= 1;
    }
  }
  while (nn != 0);


  val_elts = val.elts();
  Vec<long> tmp;
  tmp.SetLength(2*n);
  long *tmp_elts = tmp.elts();

  for (long i = 0; i < n; i++){
    tmp_elts[i] = val_elts[i]._zz_p__rep;
    tmp_elts[i + n] = 0;
  }
  reduce_t(tmp_elts);  
  for (long i = 0; i < n; i++)
    val_elts[i]._zz_p__rep = tmp_elts[i];
}

/*---------------------------------------------------------*/
/* transpose-interpolates val at n points                  */
/*---------------------------------------------------------*/
void zz_pX_Multipoint_CTFT::interpolate_t(Vec<zz_p>& f, const Vec<zz_p>& val) const{
  f.SetLength(n);
  zz_p * f_elts = f.elts();
  const zz_p * val_elts = val.elts();

  if (n == 0){
    return;
  }

  if (n == 1){
    f_elts[0] = val_elts[0];
    return;
  }


  Vec<long> wk_vec;
  wk_vec.SetLength(3*n);
  long *wk = wk_vec.elts();
  long *wk2 = wk_vec.elts() + n;
  for (long i = 0; i < n; i++)
    wk[i] = val_elts[i]._zz_p__rep;
  for (long i = 0; i < 2*n; i++)
    wk2[i] = 0;

  CRT_t(wk, wk2);
  
  Vec<zz_p> wk_zz;
  wk_zz.SetLength(n);
  zz_p * wk_zz_elts = wk_zz.elts();
  for (long i = 0; i < n; i++)
    wk_zz_elts[i]._zz_p__rep = wk[i];
  
  long nn = n;
  long k = 0;
  long aa = 1L;

  while (aa <= nn/2){
    k++;
    aa <<= 1;
  }

  do {
    interpolate_chunk_t(wk, wk_zz_elts, k, wk2);
    wk += aa;
    wk_zz_elts += aa;

    nn = nn-aa;

    aa = 1L;
    k = 0;
    while (aa <= nn/2){
      k++;
      aa <<= 1;
    }
  }
  while (nn != 0);
  
  wk = wk_vec.elts();
  for (long i = 0; i < n; i++)
    f_elts[i]._zz_p__rep = wk[i];
}
