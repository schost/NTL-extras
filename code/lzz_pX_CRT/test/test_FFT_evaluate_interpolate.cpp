#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* initializes a zz_pX_Multipoint                             */
/* check takes an extra argument, not used here               */
/*------------------------------------------------------------*/
void check(int opt){

  zz_p::FFTInit(1);
  for (long j = 98; j < 200; j+=10){
    long n = 1L << (NextPowerOfTwo(j));

    Vec<zz_p> val, pts;
    zz_pX f, X;

    zz_pX_Multipoint * ev;
    zz_pX_Multipoint_FFT evQ(n);
    ev = &evQ;

    X = 0;
    SetCoeff(X, 1, 1);
    ev->evaluate(pts, X);

    f = random_zz_pX(j);
    ev->evaluate(val, f);
    
    for (long i = 0; i < n; i++)
      if (val[i] != eval(f, pts[i]))
	cout << j << " " << i << endl;

    zz_pX g;
    ev->interpolate(g, val);
    cout << f-g << endl;
  }
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
