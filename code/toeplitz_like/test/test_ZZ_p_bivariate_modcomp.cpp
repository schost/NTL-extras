#include <NTL/lzz_pX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/vector.h>
#include <assert.h>

#include "ZZ_p_block_sylvester.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* does one test modulo p^k, with deg(p, X)=n and deg(p, Y)=m */
/* p = 2^59+1                                                 */
/*------------------------------------------------------------*/
void check(int opt, long k, long n, long m){

  ZZ p = to_ZZ(576460752303423489);
  p = power(p, k);
  ZZ_p::init(p);

  long prec = (m+1)*(n+1);

  ZZ_pX f;
  f = random_ZZ_pX(prec);

  Vec<long> type;
  type.SetLength(m+1);
  for (long i = 0; i <= m; i++)
    type[i] = n;

  Vec<ZZ_p> input, output1, output2;
  input.SetLength((m+1)*(n+1));
  for (long i = 0; i < (m+1)*(n+1); i++)
    input[i] = random_ZZ_p();

  if (opt == 1){
    ZZ_p_bivariate_modular_composition c(f, type, prec);
    output1 = c.mul_right(input);
    output2 = c.mul_right_Horners(input);
    cout << n << " " << m << " " << (output1 == output2) << endl;
  }
  else{
    cout << n << " " << m << " ";
    long NB = 1;
    double t;

    t = GetTime();
    for (long i = 0; i < NB-1; i++)
      ZZ_p_bivariate_modular_composition c(f, type, prec);
    ZZ_p_bivariate_modular_composition c(f, type, prec);
      
    cout << GetTime()-t << " ";

    t = GetTime();
    for (long i = 0; i < NB; i++)
      output1 = c.mul_right(input);
    cout << GetTime()-t << " ";

    t = GetTime();
    for (long i = 0; i < NB; i++)
      output2 = c.mul_right_Horners(input);
    cout << GetTime()-t << " ";

    cout << endl;
  }
}

void check(int opt){
  for (long i = 1; i < 200; i++)
    check(opt, 1, i, i);
}

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
