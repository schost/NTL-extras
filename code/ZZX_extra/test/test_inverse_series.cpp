#include <NTL/ZZX.h>

#include "magma_output.h"
#include "ZZX_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* check takes an extra argument, not used here               */
/*------------------------------------------------------------*/
void check(int opt){
  ZZX g, ig;
  ZZ d_g, d_ig;

  long prec = 10;
  long bits = 10;
  random(g, bits, prec);
  d_g = RandomBits_ZZ(bits-1);

  inverse_series(ig, d_ig, g, d_g, prec);
  cout << trunc(ig * g - ZZX{d_ig * d_g}, prec) << endl;
} 

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
