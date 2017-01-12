#include <NTL/lzz_pX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/vector.h>

#include "ZZ_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* multipoint evaluation over ZZ_p                            */
/* FFT points                                                 */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* bit-reversal, taken from NTL's FFT.c                       */
/*------------------------------------------------------------*/
static
long RevInc(long a, long k) {
   long j, m;

   j = k; 
   m = 1L << (k-1);

   while (j && (m & a)) {
      a ^= m;
      m >>= 1;
      j--;
   }
   if (j) a ^= m;
   return a;
}

// TODO: merge with Kevin's
static long log_order(const ZZ_p& w){
  if (w == to_ZZ_p(1))
    return 0;
  else
    return 1 + log_order(w*w);
}

/*------------------------------------------------------------*/
/* inits the array of roots of unity and the bit-rev indices  */
/*------------------------------------------------------------*/
ZZ_pX_Multipoint_FFT::ZZ_pX_Multipoint_FFT(const ZZ_p & w, long n){
  // this->k = NextPowerOfTwo(n);
  this->k = log_order(w);
    this->max_n = 1L << k;
  this->n = n;
  
    powersW.SetLength(k);
    for (long kk = 0; kk < k; kk++){
      long nn = 1L << kk;
      powersW[kk].SetLength(nn);
      ZZ_p ww = power(w, 1L << (k-1-kk));
      powersW[kk][0] = to_ZZ_p(1);
      for (long i = 1; i < nn; i++)
	powersW[kk][i] = powersW[kk][i-1] * ww;
    }

    rev.SetLength(max_n);
    long i, j;
    for (i = 0, j = 0; i < max_n; i++, j = RevInc(j, k))
      rev[i] = j;
}

/*------------------------------------------------------------*/
/* inits the array of roots of unity and the bit-rev indices  */
/*------------------------------------------------------------*/
ZZ_pX_Multipoint_FFT::ZZ_pX_Multipoint_FFT(const ZZ_p & w, const ZZ_p & c, long n) : ZZ_pX_Multipoint_FFT(w, n){
  powersC.SetLength(n);
  inv_powersC.SetLength(n);
  ZZ_p invC = 1/c;
  powersC[0] = to_ZZ_p(1);
  inv_powersC[0] = to_ZZ_p(1);
  for (long i = 1; i < n; i++){
    powersC[i] = c * powersC[i-1];
    inv_powersC[i] = invC * inv_powersC[i-1];
  }
}


/*------------------------------------------------------------*/
/* main helper function for mul_right and evaluate            */
/*------------------------------------------------------------*/
void ZZ_pX_Multipoint_FFT::evaluate_doit(Vec<ZZ_p>& val, const Vec<ZZ_p>& a) const {
  Vec<ZZ_p> dft_a;
  dft_a.SetLength(max_n);

  if (k == 0) {
    val[0] = a[0];
    return;
  }
  if (k == 1) {
    val[0] = a[0] + a[1];
    val[1] = a[0] - a[1];
    return;
  }

  // i = 0
  long maxell = 1L << (k-1);
  long s = 0;
  long t = maxell;
  for (; s < maxell; s++, t++){
    ZZ_p as = a[s];
    ZZ_p at = a[t];
    dft_a[s] = as + at;
    dft_a[t] = (as - at) * powersW[k-1][s];
  }

  // 1 <= i < k-1
  for (long i = 1; i < k-1; i++){
    maxell = (1L << (k-i-1));

    for (long j = 0; j < (1L << i); j++){
      s = j << (k-i);
      t = s + maxell;
      for (long ell = 0; ell < maxell; ell++){
        ZZ_p as = dft_a[ell + s];
        ZZ_p at = dft_a[ell + t];
        dft_a[ell + s] = as + at;
        dft_a[ell + t] = (as - at) * powersW[k-1-i][ell];
      }
    }
  }
  
  // i = k-1
  maxell = 1L << k;
  for (long s = 0; s < maxell; s += 2){
    ZZ_p as = dft_a[s];
    ZZ_p at = dft_a[s+1];
    dft_a[s] = as + at;
    dft_a[s+1] = as - at;
  }

  for (long i = 0; i < max_n; i++)
    if (rev[i] < n) // duh
      val[rev[i]] = dft_a[i];
}


/*------------------------------------------------------------*/
/* does a diagonal multiplication, then a forward FFT         */
/*------------------------------------------------------------*/
void ZZ_pX_Multipoint_FFT::evaluate(Vec<ZZ_p>& val, const ZZ_pX& f) const {

  val.SetLength(n);
  Vec<ZZ_p> a;
  a.SetLength(max_n);

  if (powersC.length() == 0)
    for (long i = 0; i < n; i++)
      a[i] = coeff(f, i);
  else
    for (long i = 0; i < n; i++)
      a[i] = coeff(f, i) * powersC[i];

  for (long i = n; i < max_n; i++)
    a[i] = to_ZZ_p(0);

  evaluate_doit(val, a);
}

/*------------------------------------------------------------*/
/* does a diagonal multiplication, then a forward FFT         */
/*------------------------------------------------------------*/
void ZZ_pX_Multipoint_FFT::mul_right(Vec<ZZ_p>& output, const Vec<ZZ_p>& input) const {

  output.SetLength(n);
  Vec<ZZ_p> a;
  a.SetLength(max_n);

  if (powersC.length() == 0)
    for (long i = 0; i < n; i++)
      a[i] = input[i];
  else
    for (long i = 0; i < n; i++)
      a[i] = input[i] * powersC[i];

  for (long i = n; i < max_n; i++)
    a[i] = to_ZZ_p(0);

  evaluate_doit(output, a);
}

/*------------------------------------------------------------*/
/* does a forward FFT, then a diagonal multiplication         */
/*------------------------------------------------------------*/
void ZZ_pX_Multipoint_FFT::mul_left(Vec<ZZ_p>& output, const Vec<ZZ_p>& input) const {

  output.SetLength(n);
  Vec<ZZ_p> a;
  a.SetLength(max_n);

  for (long i = 0; i < n; i++)
    a[i] = input[i];
  for (long i = n; i < max_n; i++)
    a[i] = to_ZZ_p(0);

  evaluate_doit(output, a);

  if (powersC.length() != 0)
    for (long i = 0; i < n; i++)
      output[i] = output[i] * powersC[i];
}


/*------------------------------------------------------------*/
/* finds a root of unity of order s mod p                     */
/* assumes s is a power of 2                                  */
/*------------------------------------------------------------*/
long find_root_of_unity(long p, long s){
  long ss = s;
  long log = 0;
  while ((ss & 1) == 0){
    log++;
    ss = ss >> 1;
  }
  if (ss != 1)
    LogicError("FFT init with non power of 2");

  zz_pPush push;
  zz_p::UserFFTInit(p);
  long w;
  IsFFTPrime(p, w);
  zz_p zww = to_zz_p(w);

  long order = 0;
  while (zww != to_zz_p(1)){
    zww *= zww;
    order++;
  }

  zww = power(to_zz_p(w), 1L << (order-log));

  return rep(zww);
}

/*------------------------------------------------------------*/
/* suppose that omega is a root of unity of order s mod p     */  
/* lifts omega to a root of unity Omega mod p^k               */
/*------------------------------------------------------------*/
void lift_root_of_unity(ZZ& Omega, long omega, long s, long p, long k){
  ZZ_pPush push;
  Omega = to_ZZ(omega);

  if (k == 1)
    return;

  ZZ pk = to_ZZ(p);

  long prec = 1;
  do{
    pk = pk*pk;
    ZZ_p::init(pk);
    ZZ_p s_mod = to_ZZ_p(s);
    ZZ_p Omega_mod = to_ZZ_p(Omega);
    ZZ_p lift = Omega_mod*(to_ZZ_p(1) - (power(Omega_mod, s)-1)/s_mod);
    Omega = rep(lift);
    prec = 2*prec;
  } while (prec < k);

  Omega = Omega % power(to_ZZ(p), k);
}
