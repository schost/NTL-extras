#include <NTL/lzz_pX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/vector.h>
#include <assert.h>

#include "ZZ_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* initializes a ZZ_pX_Multipoint                             */
/* check takes an extra argument, not used here               */
/*------------------------------------------------------------*/
void check(int opt){

  long p, w, s;
  s = 1L << 1;
  zz_p::FFTInit(1);
  p = zz_p::modulus();
  w = find_root_of_unity(p, s);

  long k = 100;
  ZZ Omega;
  lift_root_of_unity(Omega, w, s, p, k);

  ZZ_p::init(power(to_ZZ(p), k));
  ZZ_p W = to_ZZ_p(Omega);
  assert (power(W, s >> 1) != to_ZZ_p(1));
  assert (power(W, s) == to_ZZ_p(1));
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
