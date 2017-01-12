#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_CRT.h"
#include "magma_output.h"
#include "lzz_pXY.h"
#include "lzz_pEX_augmented.h"

NTL_CLIENT

void local_crt(zz_pX & res, const zz_pX & a0, const zz_pX& a1, const zz_pX& a2, zz_pX_CRT& crt){
  Vec<zz_pX> rems;
  rems.SetLength(3);
  rems[0] = a0;
  rems[1] = a1;
  rems[2] = a2;
  crt.combine(res, rems);
}

/*------------------------------------------------------------*/
/* check takes an extra argument, not used here               */
/*------------------------------------------------------------*/
void check(int opt){

  long p = 1125899906842679;
  zz_p::init(p);

  zz_pX f0, f1, f2;
  long d = 5;

  f0 = random_zz_pX(d);
  SetCoeff(f0, d, 1);
  f1 = random_zz_pX(d);
  SetCoeff(f1, d, 1);
  f2 = random_zz_pX(d);
  SetCoeff(f2, d, 1);

  Vec<zz_pX> mods;
  mods.SetLength(3);
  mods[0] = f0;
  mods[1] = f1;
  mods[2] = f2;
  zz_pX_CRT crt(mods);

  Vec<zz_pX> coeffs;
  coeffs.SetLength(4);
  local_crt(coeffs[0], to_zz_pX(0), random_zz_pX(d), random_zz_pX(d), crt);
  local_crt(coeffs[1], to_zz_pX(0), random_zz_pX(d), random_zz_pX(d), crt);
  local_crt(coeffs[2], to_zz_pX(0), to_zz_pX(0), random_zz_pX(d), crt);
  local_crt(coeffs[3], to_zz_pX(0), to_zz_pX(0), random_zz_pX(d), crt);

  zz_pEX_augmented T(crt.master(), coeffs);
  Vec<zz_pEX_augmented> Tsplit;
  make_monic_and_split(Tsplit, T);

  magma_init();
  magma_init_bi();
  magma_assign(T, "T");
  magma_assign(Tsplit[0], "T0");
  magma_assign(Tsplit[1], "T1");
  magma_assign(Tsplit[2], "T2");
  
  cout << "print Ideal(T) eq Ideal(T0) * Ideal(T1) * Ideal(T2);\n";
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
