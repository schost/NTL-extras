#include <NTL/ZZX.h>
#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_CRT.h"
#include "magma_output.h"
#include "ZZXY.h"
#include "lzz_pXY.h"
#include "lzz_pEX_augmented.h"

NTL_CLIENT

void check(const string & name){

  long p = 1125899906842679;
  zz_p::init(p);

  ZZXY FZ = ZZXY();
  FZ.read_from_file(name);

  zz_pXY F, Fy;
  conv(F, FZ);
  diffY(Fy, F);
  
  magma_init();
  magma_init_bi();
  magma_assign_bi(F, "F");
  magma_assign_bi(Fy, "Fy");
  
  cout << "I:=Ideal([F,Fy]);\n";
  
  Vec<zz_pEX_augmented> sols;
  solve(sols, F, Fy);

  cout << "sol:=[];\n";
  for (long i = 0; i < sols.length(); i++)
    magma_assign(sols[i], "sol[" + to_string(i) + "+1]");
  
  if (sols.length() != 0)
    cout << "print {GroebnerBasis(y) :y in &cat [ProbableRadicalDecomposition(Ideal(x)): x in sol]} eq {GroebnerBasis(y) :y in ProbableRadicalDecomposition(Ideal(I))};\n";
  else 
    cout << "print Dimension(I) ne 0;\n";
  
}  

int main(int argc, char ** argv){

  string name;
  if (argc > 1)
    name = argv[1];
  else
    name = "../data/test_OL_16_009.txt.ntl";

  check(name);
  return 0;
}
