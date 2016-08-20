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

  zz_pX f0, g0, f1, g1;

  long d = 100;

  g0 = random_zz_pX(d);
  f0 = g0+1;
  g1 = random_zz_pX(d/2);
  f1 = 0;
  
  Vec<zz_pX> vecF, vecG;
  vecF.SetLength(2);
  vecF[0] = f0;
  vecF[1] = f1;
  vecG.SetLength(2);
  vecG[0] = g0;
  vecG[1] = g1;
  zz_pX_CRT crt(vecG);
  
  zz_pX f, g; 
  g = g0*g1;
  crt.combine(f, vecF);
  cout << (f - f0) % g0 << " ";
  cout << (f - f1) % g1 << endl;

  zz_pX zero, unit, inv;
  split(zero, unit, inv, f, g);

  cout << zero - g1/LeadCoeff(g1) << " ";
  cout << unit - g0/LeadCoeff(g0) << " ";
  cout << inv*f % unit << endl;
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
