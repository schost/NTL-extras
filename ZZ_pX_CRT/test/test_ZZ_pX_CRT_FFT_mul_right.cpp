#include <NTL/lzz_pX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/vector.h>
#include <assert.h>

#include "ZZ_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* initializes a ZZ_pX_Multipoint                             */
/* check takes an extra argument, not used here               */
/*------------------------------------------------------------*/
void check(int opt){

  zz_p::FFTInit(1);
  long p = zz_p::modulus();
  long k = 11;
  ZZ pZZ = to_ZZ(p);
  ZZ pk = power(pZZ, k);
  ZZ_p::init(pk);

  for (long i = 1; i < 10000; i = 2*i+3){
    long order = 1L << NextPowerOfTwo(i);
    long w = find_root_of_unity(p, order);
    ZZ W;
    lift_root_of_unity(W, w, order, p, k);
    ZZ_p omega = to_ZZ_p(W);
    ZZ_p c = random_ZZ_p();

    ZZ_pX_Multipoint_FFT fft(omega, c, i);
    ZZ_pX f = random_ZZ_pX(i);
    Vec<ZZ_p> coeffs;
    coeffs.SetLength(i);
    for (long k = 0; k < i; k++)
      coeffs[k] = coeff(f, k);

    Vec<ZZ_p> val;
    
    if (opt == 1){
      fft.mul_right(val, coeffs);
      Vec<ZZ_p> val2;
      val2.SetLength(i);
      for (long j = 0; j < i; j++)
	val2[j] = eval(f, c*power(omega, j));
      assert (val == val2);
    }
    else{
      double t;
      cout << i << " ";
      
      t = GetTime();
      for (long j = 0; j < 100; j++)
	fft.evaluate(val, f);
      t = GetTime() - t;
      cout << t << " ";
      
      t = GetTime();
      ZZ_pX g = random_ZZ_pX(i);
      for (long j = 0; j < 100; j++)
	ZZ_pX h = f*g;
      t = GetTime() - t;
      cout << t << " ";
      
      cout << endl;
    }
  }
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
