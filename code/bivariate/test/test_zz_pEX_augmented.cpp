#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_CRT.h"
#include "magma_output.h"
#include "lzz_pXY.h"
#include "lzz_pEX_augmented.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* check takes an extra argument, not used here               */
/*------------------------------------------------------------*/
void check(int opt){

  long p = 1125899906842679;
  zz_p::init(p);

  zz_pX t1;
  zz_pXY t2;

  long d1 = 10, d2 = 3;

  t1 = random_zz_pX(d1);
  SetCoeff(t1, d1, 1);

  random_zz_pXY(t2, d1, d2);

  zz_pEX_augmented T(t1, t2);
  magma_init();
  magma_init_bi();
  magma_assign(T, "t");
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
