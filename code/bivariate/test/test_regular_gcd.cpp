#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_CRT.h"
#include "magma_output.h"
#include "ZZXY.h"
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

  Vec<pair_zz_pEX_augmented> T;
  T.SetLength(9);

  random_pair_monic_zz_pEX_augmented(T[0], 5, -1, -1);
  random_pair_monic_zz_pEX_augmented(T[1], 5, -1, 4);
  random_pair_monic_zz_pEX_augmented(T[2], 5, -1, 8);

  random_pair_monic_zz_pEX_augmented(T[3], 5, 3, 1);
  random_pair_monic_zz_pEX_augmented(T[4], 5, 3, 4);
  random_pair_monic_zz_pEX_augmented(T[5], 5, 3, 8);

  random_pair_monic_zz_pEX_augmented(T[6], 5, 6, 1);
  random_pair_monic_zz_pEX_augmented(T[7], 5, 6, 4);
  random_pair_monic_zz_pEX_augmented(T[8], 5, 6, 8);

  for (long i = 0; i < 9; i++){
    zz_pEPush push(T[i].T1);
    zz_pEX tmp = random_zz_pEX(3+i);
    T[i].T20 *= tmp * random_zz_pEX(3); 
    T[i].T21 *= tmp * random_zz_pEX(3); 
  }

  pair_zz_pEX_augmented Tall;
  combine(Tall, T);

  Vec<zz_pEX_augmented> Tsplit;
  regular_gcd(Tsplit, Tall);

  for (long i = 0; i < Tsplit.length(); i++)
    cout << deg(Tsplit[i].T1X) << " " << deg(Tsplit[i].T2) << endl;
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
