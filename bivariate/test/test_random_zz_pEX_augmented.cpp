#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_CRT.h"
#include "magma_output.h"
#include "lzz_pXY.h"
#include "lzz_pEX_augmented.h"

NTL_CLIENT

zz_pX random_monic_zz_pX(long d){
  zz_pX f = random_zz_pX(d);
  SetCoeff(f, d, 1);
  return f;
}

/*------------------------------------------------------------*/
/* check takes an extra argument, not used here               */
/*------------------------------------------------------------*/
void check(int opt){

  long p = 1125899906842679;
  zz_p::init(p);

  zz_pEX_augmented T;
  random_monic_zz_pEX_augmented(T, 3, 5);
  
  magma_init();
  magma_init_bi();
  magma_assign(T, "T");
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
