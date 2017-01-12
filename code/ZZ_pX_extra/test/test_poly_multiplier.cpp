#include <NTL/ZZ_pX.h>
#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "ZZ_pX_extra.h"
#include "magma_output.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* does a multipoint evaluation                               */
/*------------------------------------------------------------*/
void check(int opt){
  ZZ p = to_ZZ("234898888888888888888888888888888888888888881111111111");
  ZZ_p::init(p);

  long NB = 1000;
  
  for (long i = 1; i < 1000; i++){
    long dA = i;
    long dB = i;
    ZZ_pX A, B, C, C2;
    A = random_ZZ_pX(dA);
    B = random_ZZ_pX(dB);
    
    cout << i << " ";
    double t;
    
    t = GetTime();
    ZZ_pX_poly_multiplier multA(A, dB);
    cout << GetTime()-t << " ";

    t = GetTime();
    for (long j = 0; j < NB; j++)
      multA.mul(C, B);
    cout << GetTime()-t << " ";

    t = GetTime();
    for (long j = 0; j < NB; j++)
      C2 = A*B;
    cout << GetTime()-t << " ";

    cout << (C-C2) << endl;
  }
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
