#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* computes the FFT tables                                    */
/*------------------------------------------------------------*/
void check(int opt){
  zz_p::FFTInit(1);
  zz_pX_Multipoint_CTFT mult1;
  mult1.init_multipliers(5);

  zz_p::FFTInit(0);
  zz_pX_Multipoint_CTFT mult2;
  mult1.init_multipliers(5);
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
