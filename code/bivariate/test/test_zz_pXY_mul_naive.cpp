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

  zz_pXY F, G, H, I, J;

  long d = 3;

  random_zz_pXY(F, d, d);
  mul_naive(H, G, F);
  random_zz_pXY(I, 2*d, d);
  mul_naive(J, I, F);

  magma_init();
  magma_init_bi();
  magma_assign_bi(F, "F");
  magma_assign_bi(G, "G");
  magma_assign_bi(H, "H");
  magma_assign_bi(I, "I");
  magma_assign_bi(J, "J");
  cout << "print H eq F*G;\n";
  cout << "print J eq F*I;\n";
  
  
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
