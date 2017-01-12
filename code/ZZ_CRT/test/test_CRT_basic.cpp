#include <NTL/ZZ.h>
#include <NTL/vector.h>

#include "magma_output.h"
#include "ZZ_CRT.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/* does a multimod                                            */
/* check takes an extra argument, not used here               */
/*------------------------------------------------------------*/
void check(int opt){
  Vec<long> moduli;

  for (long j = 10; j < 11; ){
    long base = 1L << 59;
    moduli.SetLength(j);
    for (long i = 0; i < moduli.length(); i++){
      base = NextPrime(base+1);
      moduli[i] = base;
    }
    
    ZZ_CRT_crt_basic crt(moduli);
    ZZ_CRT_rem_basic rem(moduli);

    Vec<long> rems;
    ZZ a = to_ZZ("1725436586697640946858688965569256363112777243042596638790631055949824");
    for (long i = 0; i < moduli.length(); i++)
      a *= (3*to_ZZ(moduli[i])/2);
    a = -a;

    rem.eval(rems, a);
    ZZ b;
    crt.crt(b, rems);

    if (0){
      cout << j << endl;
    }
    else{
      cout << "a:=" << a << ";\n";
      magma_assign(moduli, "m");
      magma_assign(rems, "r");
      cout << "b:=" << b << ";\n";
      cout << "print b eq (a mod &*m);\n";
    }
   } 
  
} 

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
