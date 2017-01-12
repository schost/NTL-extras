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

  zz_pXY F, G, I, J;

  long d = 3;

  random_zz_pXY(F, d, d);
  random_zz_pXY(G, 2*d, 2*d);
  random_zz_pXY(I, 2*d, 0);
  mul_naive(G, F, G);
  mul_naive(G, G, I);
  mul_naive(F, F, I);
  
  GCD(J, F, G);

  magma_init();
  magma_init_bi();
  magma_assign_bi(F, "F");
  magma_assign_bi(G, "G");
  magma_assign_bi(J, "J");
  cout << "print Degree(GCD(F, G) div J) eq 0;\n";
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
