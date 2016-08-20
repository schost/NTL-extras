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

  for (long d = 0; d < 10; d++){
    Vec<zz_pX> F;
    F.SetLength(d);
    for (long i = 0; i < d; i++){
      F[i] = random_zz_pX(d);
    }
    
    zz_pXY Fpoly(F);
    zz_pXY dFpoly;
    diffY(dFpoly, Fpoly);

    magma_init();
    magma_init_bi();
    magma_assign_bi(F, "F");
    magma_assign_bi(dFpoly, "Fy");
    cout << "print Fy eq Derivative(F, Y);\n";
  }
    
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
