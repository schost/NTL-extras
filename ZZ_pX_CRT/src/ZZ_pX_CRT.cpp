#include <NTL/ZZ_pX.h>
#include <NTL/vector.h>

#include "ZZ_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* multipoint evaluation over zz_p                            */
/* FFT points                                                 */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/


ZZ_pX_Multipoint_FFT::ZZ_pX_Multipoint_FFT(const ZZ_p & w, long n){
    this->k = NextPowerOfTwo(n);
    this->s = 1L << k;
    this->n = n;
    if (ZZ_pInfo->p_info == NULL){
      cerr << "Attempt to init a ZZ_pX_Multipoint_FFT without ZZ_p::FFTInit\n";
      exit(-1);
    }

    powersW.SetLength(k+1);
    for (long kk = 0; kk < k; kk++){
      long nn = 1L << k;
      powersW[k].SetLength(nn);
      ZZ_p ww = power(w, 1L << (k-kk));
      powersW[kk][0] = to_ZZ_p(1);
      for (i = 1; i < nn; i++)
	powersW[k][i] = powersW[k][i-1] * ww;
    }
}

/*------------------------------------------------------------*/
/* does a forward FFT                                         */
/*------------------------------------------------------------*/
void ZZ_pX_Multipoint_FFT::evaluate(Vec<ZZ_p>& val, const ZZ_pX& f) const {

  val.SetLength(n);
  Vec<ZZ_p> a, dft_a;
  a.SetLength(s);
  dft_a.SetLength(s);
  for (long i = 0; i < n; i++)
    a[i] = coeff(f, i);
  for (long i = n; i < s; i++)
    a[i] = to_ZZ_p(0);

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
    as = a[s];
    at = a[t];
    dft_a[s] = as + at;
    dft_a[t] = (as - at) * w[k][s];
  }
  // 1 <= i < k-1
  for (i = 1; i < k-1; i++){

    for (j = 0; j < (1L << i); j++){
      maxell = (1L << (k-i-1));
      s = j << (k-i);
      t = s + maxell;

      for (ell = 0; ell < maxell; ell++){
        as = dft_a[ell + s];
        at = dft_a[ell + t];
        dft_a[ell + s] = as + at;
        dft_a[ell + t] = (as - at) * wi[k-i][ell];
      }
    }
  }
  
  // i = k-1
  maxell = 1L << k;
  for (s = 0; s < maxell; s += 2){
    as = dft_a[s];
    at = dft_a[s+1];
    dft_a[s] = as + at;
    dft_a[s+1] = as - at;
  }

  for (long i = 0; i < n; i++)
    val[i] = dft_a[i];

}
