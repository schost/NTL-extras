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

  zz_pX t1, t1a, t1b, t1c;
  zz_pXY t20, t21;

  long d1 = 5, d2 = 3;

  t1a = random_zz_pX(d1);
  SetCoeff(t1a, d1, 1);
  t1b = random_zz_pX(d1);
  SetCoeff(t1b, d1, 1);
  t1c = random_zz_pX(d1);
  SetCoeff(t1c, d1, 1);
  t1 = t1a*t1b*t1c;

  random_zz_pXY(t20, deg(t1), d2);
  random_zz_pXY(t21, deg(t1), d2+3);
  pair_zz_pEX_augmented T(t1, t20, t21);
  magma_init();
  magma_init_bi();
  magma_assign(T, "t");

  Vec<zz_pX> mods;
  mods.SetLength(3);
  mods[0] = t1a;
  mods[1] = t1b;
  mods[2] = t1c;
  zz_pX_CRT crt(mods);

  Vec<pair_zz_pEX_augmented> Tmod;
  multimod(Tmod, T, mods);
  
  magma_assign(Tmod[0], "tmod0");
  magma_assign(Tmod[1], "tmod1");
  magma_assign(Tmod[2], "tmod2");

  magma_assign(t1a, "X", "u0");
  magma_assign(t1b, "X", "u1");
  magma_assign(t1c, "X", "u2");

  cout << "print (Ideal([t[1],t[2]]) + Ideal([u0])) eq Ideal([tmod0[1],tmod0[2]]);\n";
  cout << "print (Ideal([t[1],t[3]]) + Ideal([u0])) eq Ideal([tmod0[1],tmod0[3]]);\n";
  cout << "print (Ideal([t[1],t[2]]) + Ideal([u1])) eq Ideal([tmod1[1],tmod1[2]]);\n";
  cout << "print (Ideal([t[1],t[3]]) + Ideal([u1])) eq Ideal([tmod1[1],tmod1[3]]);\n";
  cout << "print (Ideal([t[1],t[2]]) + Ideal([u2])) eq Ideal([tmod2[1],tmod2[2]]);\n";
  cout << "print (Ideal([t[1],t[3]]) + Ideal([u2])) eq Ideal([tmod2[1],tmod2[3]]);\n";

 
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
