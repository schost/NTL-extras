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
  long d = 7;
  long b = 2;

  ZZXY F;
  random(F, b, d, d);

  ZZX F0;
  random(F0, b, F.degY()-1);
  ZZX factor;
  ZZ root{22};
  SetCoeff(factor, 1, 1);
  SetCoeff(factor, 0, -root);
  F0 = F0 * factor;

  for (long i = 0; i <= F.degY(); i++)
    F.coeffX[i] = (F.coeffX[i] << 1) + coeff(F0, i);
    
  
  ZZX g, h;
  ZZ den_g, den_h;

  F.series_solution(g, den_g, root, d*d);
  F.eval(h, den_h, g, den_g, d*d);
  cout << h << endl;
  
  // //F.eval(h, g, d*d);

  // magma_init_bi_QQ();
  // magma_init_QQX();
  // magma_assign(F, "F");
  // magma_assign(g, "XX", "g");
  // cout << "gr:=1/" << den_g << "*g;\n";
  // magma_assign(h, "XX", "h");
  // cout << "hr:=1/" << den_h << "*h;\n";
  // cout << "(hr-Evaluate(F, [gr,XX])) mod XX^" << d*d << ";\n";
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
