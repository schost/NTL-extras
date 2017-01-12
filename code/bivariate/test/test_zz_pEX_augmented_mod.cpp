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

  zz_pX t1, t1a, t1b;
  zz_pXY t2;

  long d1 = 5, d2 = 3;

  t1a = random_zz_pX(d1);
  SetCoeff(t1a, d1, 1);
  t1b = random_zz_pX(d1);
  SetCoeff(t1b, d1, 1);
  t1 = t1a*t1b;

  random_zz_pXY(t2, deg(t1), d2);
  zz_pEX_augmented T(t1, t2);
  magma_init();
  magma_init_bi();
  magma_assign(T, "t");

  zz_pEX_augmented Tmod;
  mod(Tmod, T, t1a);
  magma_assign(Tmod, "tmod");

  magma_assign(t1a, "X", "u");
  
  cout << "print (Ideal(t) + Ideal([u])) eq Ideal(tmod);\n";
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
