#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_CRT.h"
#include "magma_output.h"
#include "ZZXY.h"
#include "lzz_pXY.h"
#include "lzz_pEX_augmented.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* check takes an extra argument, not used here               */
/*------------------------------------------------------------*/
void check(int opt){

  long p = 1125899906842679;
  zz_p::init(p);

  zz_pXY F, G;
  random_zz_pXY(F, 10, 10);
  random_zz_pXY(G, 10, 10);
  magma_init();
  magma_init_X();
  magma_init_Y();
  magma_assign(F, "F");
  magma_assign(G, "G");
  
  zz_pX res;
  resultant(res, F, G);
  magma_assign(res, "res");
  cout << "print res eq Resultant(F,G);\n";
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
