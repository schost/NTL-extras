#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "magma_output.h"
#include "lzz_pXY.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* initializes output, print some poly                        */
/* check takes an extra argument, not used here               */
/*------------------------------------------------------------*/
void check(int opt){

  long p = 1125899906842679;
  zz_p::init(p);

  for (long d = 10; d < 110; d+=10){

    Vec<zz_pX> F;
    F.SetLength(d);
    for (long i = 0; i < d; i++){
      F[i] = random_zz_pX(d);
    }
    
    cout << d << endl;

    zz_pXY Fpoly(F);

    zz_pX R = random_zz_pX(d*d);
    SetCoeff(R, d*d, 1);
    zz_pX S = random_zz_pX(d*d); 

    double t;

    t = GetTime();
    zz_pX T;
    reduce(T, R, S, F);
     cout << GetTime()-t << " ";

    t = GetTime();
    zz_pX T_naive;
    reduce_naive(T_naive, R, S, F);
    cout << GetTime()-t << " ";

    cout << endl;

    if (T != T_naive)
      cout << "print false;\n";

    // magma_init();
    // magma_init_bi();
    // magma_init_X();

    // magma_assign_bi(F, "F");
    // magma_assign(R, "R");
    // magma_assign(S, "S");
    // magma_assign(T, "T");

    // cout << "Q<_x>:=quo<U|R>;\n";
    // cout << "_y:=Q!S;\n";
    // cout << "rem:=Evaluate(F, [_y,_x]);\n";
    // cout << "rem2:=Q!T;\n";
    // cout << "print rem eq rem2;\n";
  }
    
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
