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


  ZZXY FZ = ZZXY();
  FZ.read_from_file(name);
  ZZXY FZy;
  diffY(FZy, FZ);
  
  magma_init_bi_QQ();
  magma_assign(FZ, "F");
  magma_assign(FZy, "Fy");

  Vec<ZZ_bivariate_regular_chain> sols;
  solve(sols, FZ, FZy);

  magma_assign(sols, "sols");

  cout << "p:=NextPrime(2^100);\n";
  cout << "Mp:=PolynomialRing(GF(p),2);\n";
  cout << "f:=Evaluate(F,[Mp.1, Mp.2]);\n";
  cout << "g:=Evaluate(Fy,[Mp.1, Mp.2]);\n";
  cout << "sols_p:=[ [ Evaluate(popo,[Mp.1, Mp.2]) : popo in ff] : ff in sols];\n";
  cout << "S:=0;\n";
  cout << "for I in sols_p do\n";
  cout << "    if not IsSquarefree(I[1]) then print \"output not radical\"; end if;\n";
  cout << "    for J in sols_p do\n";
  cout << "	if I ne J and Degree(GCD(I[1], J[1])) ge 1 then print \"output not coprime\"; end if;\n";
  cout << "    end for;\n";
  cout << "    U1:=UnivariatePolynomial(I[1]);\n";
  cout << "    Q:=quo<Parent(U1) | U1>;\n";
  cout << "    den:=1/Q!Derivative(U1);\n";
  cout << "    new:=&+[ Evaluate( Parent(U1)!(den*Q!UnivariatePolynomial(Coefficient(I[2], Mp.1, i))), Mp.2 )*Mp.1^i : i in [0..Degree(I[2], Mp.1)]];\n";
  cout << "    Q:=quo<Mp|[I[1], new]>;\n";
  cout << "    if Q!f ne 0 or Q!g ne 0 then print \"polynomials do not reduce to zero\"; end if;\n";
  cout << "    S:=S+Dimension(Q);\n";
  cout << "end for;\n";
  cout << "Mp:=PolynomialRing(GF(p),2,\"grevlex\");\n";
  cout << "f:=Evaluate(F,[Mp.1, Mp.2]);\n";
  cout << "g:=Evaluate(Fy,[Mp.1, Mp.2]);\n";
  cout << "R:=ProbableRadicalDecomposition(Ideal([f,g]));\n";
  cout << "T:=&+[Dimension(quo<Mp|r>) : r in R];\n";
  cout << "if S ne T then print \"missing solutions\"; end if;\n";
}  

int main(int argc, char ** argv){

  string name;
  if (argc > 1)
    name = argv[1];
  else
    name = "../data/test_OL_16_009.txt.ntl";

  double t = GetTime();
  check(name);
  cerr << "-----time: " << GetTime()-t << endl;
  return 0;
}
