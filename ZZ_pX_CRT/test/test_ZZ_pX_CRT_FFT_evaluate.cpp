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
  long k = 1;
  ZZ pZZ = to_ZZ(p);
  ZZ pk = power(pZZ, k);
  ZZ_p::init(pk);
  

  for (long i = 8; i < 10; i = 2*i){
    long w = find_root_of_unity(p, i);
    ZZ W;
    lift_root_of_unity(W, w, i, p, k);
    ZZ_p omega = to_ZZ_p(W);
   

    ZZ_pX_Multipoint_FFT fft(omega, i);
    ZZ_pX f = random_ZZ_pX(i);
    Vec<ZZ_p> val;
    
    if (opt == 1){
      fft.evaluate(val, f);
      Vec<ZZ_p> val2;
      val2.SetLength(i);
      for (long j = 0; j < i; j++)
	val2[j] = eval(f, power(omega, j));
      cout << val << endl;
      cout << val2 << endl;
      cout << endl;
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
