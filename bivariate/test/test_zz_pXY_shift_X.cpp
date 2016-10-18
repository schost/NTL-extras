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

  zz_pXY F, G;

  long d = 3;

  zz_p c = random_zz_p();
  random_zz_pXY(F, 2*d, d);
  shift_X(G, F, c);


  magma_init();
  magma_init_bi();
  magma_assign_bi(F, "F");
  magma_assign_bi(G, "G");
  cout << "c:=k!" << c << ";\n";
  cout << "print G eq Evaluate(F, [Y,X+c]);\n";
  
  
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
