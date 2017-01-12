#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "magma_output.h"
#include "ZZX_extra.h"
#include "ZZXY.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* evaluate F(x,g(x)) mod x^(d^2)                             */
/* check takes an extra argument, not used here               */
/*------------------------------------------------------------*/
void check(int opt){
  long d = 10;
  long b = 10;

  ZZXY F;
  random(F, 10, d, d);
  ZZX g, h;
  ZZ den_g, den_h;
  random(g, b, d*d);
  den_g = RandomBits_ZZ(b-1);

  F.eval(h, den_h, g, den_g, d*d);
  //F.eval(h, g, d*d);

  magma_init_bi_QQ();
  magma_init_QQX();
  magma_assign(F, "F");
  magma_assign(g, "XX", "g");
  cout << "gr:=1/" << den_g << "*g;\n";
  magma_assign(h, "XX", "h");
  cout << "hr:=1/" << den_h << "*h;\n";
  cout << "(hr-Evaluate(F, [gr,XX])) mod XX^" << d*d << ";\n";
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
