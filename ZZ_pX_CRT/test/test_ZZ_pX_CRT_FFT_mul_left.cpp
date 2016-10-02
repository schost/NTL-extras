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
    Vec<ZZ_p> input, output;
    input.SetLength(i);
    for (long k = 0; k < i; k++)
      input[k] = random_ZZ_p();
    Vec<ZZ_p> val;
    
    fft.mul_left(val, input);
    Vec<ZZ_p> val2;
    val2.SetLength(i);
    for (long j = 0; j < i; j++){
      val2[j] = 0;
      for (long k = 0; k < i; k++)
	val2[j] = val2[j] + input[k] * power(c*power(omega, k), j);
    }
    assert (val == val2);

  }
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
