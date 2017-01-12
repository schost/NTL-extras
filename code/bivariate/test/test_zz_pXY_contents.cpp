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

  zz_pXY F, I, J;

  long d = 3;

  random_zz_pXY(F, d, d);
  random_zz_pXY(I, 2*d, 0);
  mul_naive(J, I, F);

  zz_pX c;
  contents(c, J);
  
  magma_init();
  magma_init_bi();
  magma_assign_bi(J, "J");
  magma_assign(c, "X", "c");
   
  cout << "print c eq GCD(Coefficients(J, Y));\n";
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
