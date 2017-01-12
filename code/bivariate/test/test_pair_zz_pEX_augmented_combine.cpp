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

  pair_zz_pEX_augmented Ta;
  random_pair_monic_zz_pEX_augmented(Ta, 3, -1, 2);

  pair_zz_pEX_augmented Tb;
  random_pair_monic_zz_pEX_augmented(Tb, 4, 0, 1);
  
  magma_init();
  magma_init_bi();
  magma_assign(Ta, "Ta");
  magma_assign(Tb, "Tb");

  Vec<pair_zz_pEX_augmented> lT;
  lT.SetLength(2);
  lT[0] = Ta;
  lT[1] = Tb;

  pair_zz_pEX_augmented T;
  combine(T, lT);

  magma_assign(T, "T");
  cout << "print M!(quo<M|[Ta[1]]>!T[2]) eq Ta[2];\n";
  cout << "print M!(quo<M|[Tb[1]]>!T[2]) eq Tb[2];\n";
  cout << "print M!(quo<M|[Ta[1]]>!T[3]) eq Ta[3];\n";
  cout << "print M!(quo<M|[Tb[1]]>!T[3]) eq Tb[3];\n";

}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
