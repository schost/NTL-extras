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

  zz_pXY Fpoly;

  long d = 3;

  random_zz_pXY(Fpoly, d, d);
  magma_init();
  magma_init_X();
  magma_init_Y();
  magma_assign(Fpoly, "F");
  magma_init_bi();
  magma_assign_bi(Fpoly, "FF");
  cout << "print Evaluate(FF,[y,x])-F;\n";
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
