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

  for (long j = 1; j < 1000; j++){
    long base = 1L << 59;
    moduli.SetLength(j);
    for (long i = 0; i < moduli.length(); i++){
      base = NextPrime(base+1);
      moduli[i] = base;
    }
    
    ZZ_CRT crt(moduli);


    ZZ a = to_ZZ("1725436586697640946858688965569256363112777243042596638790631055949824");
    for (long i = 0; i < moduli.length(); i++)
      a *= (3*to_ZZ(moduli[i])/2);

    ZZX apol;
    for (long i = 0; i < 10; i++)
      SetCoeff(apol, i, RandomBnd(a));
      
    Vec<zz_pX> rem;
    crt.multimod(rem, apol);

    if (1)
      cout << j << endl;
    else{
      magma_init_ZZX();
      magma_assign(apol, "a");
      magma_assign(moduli, "m");
      magma_assign(rem, "XX", "r");
      cout << "print " << j << ";\n";
      cout << "print &and [PolynomialRing(Integers(m[i]))!(r[i]-a) eq 0 : i in [1..#m]];\n";
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
