#include <NTL/ZZ_pX.h>
#include <NTL/lzz_pX.h>

NTL_CLIENT

void to_modular_rep(vec_long& x, const ZZ_p& a, const ZZ_pFFTInfoT *FFTInfo, ZZ_pTmpSpaceT *TmpSpace){
   FFTInfo->rem_struct.eval(&x[0], rep(a),  TmpSpace->rem_tmp_vec);
}

// NOTE: a gets destroyed
void from_modular_rep(ZZ_p& x, Vec<long>& avec, const ZZ_pFFTInfoT *FFTInfo, ZZ_pTmpSpaceT *TmpSpace){
   NTL_ZZRegister(t);
   long * NTL_RESTRICT a = avec.elts();

   if (FFTInfo->crt_struct.special()) {
       FFTInfo->crt_struct.eval(t, a, TmpSpace->crt_tmp_vec);
      x.LoopHole() = t;
      return;
   }

   long nprimes = FFTInfo->NumPrimes;
   const long *u = FFTInfo->u.elts();
   const long *prime = FFTInfo->prime.elts();
   const mulmod_precon_t  *uqinv = FFTInfo->uqinv.elts();
   const double *prime_recip = FFTInfo->prime_recip.elts();
      
   double y = 0.0;

   for (long i = 0; i < nprimes; i++) {
      long r = MulModPrecon(a[i], u[i], prime[i], uqinv[i]);
      a[i] = r;
      y += double(r)*prime_recip[i];
   }

   long q = long(y + 0.5);

   FFTInfo->crt_struct.eval(t, a, TmpSpace->crt_tmp_vec);

   MulAddTo(t, FFTInfo->MinusMModP, q);
   FFTInfo->reduce_struct.eval(x.LoopHole(), t);
}
