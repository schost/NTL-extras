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

  zz_pXY F;
  long d = 10;
  random_zz_pXY(F, 2*d, d);
  zz_p c = random_zz_p();
  zz_pX Fx;
  evaluate(Fx, F, c);

  magma_init();
  magma_init_bi();
  magma_assign_bi(F, "F");
  magma_assign(Fx, "Y", "Fx");
  cout << "c:=" << c << ";\n";
  cout << "print Fx eq Evaluate(F, [Y,c]);\n";
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
