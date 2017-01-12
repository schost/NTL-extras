#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>

#include "magma_output.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* initializes a zz_pX_Multipoint                             */
/* check takes an extra argument, not used here               */
/*------------------------------------------------------------*/
void check(int opt){

  long i = 30;
  long j = 500;
  Mat<zz_pX> a, b, c1, c2;

  cout << i << " " << j << " ";

  double t;
  zz_p::init(1125899906842679);

  random_mat_zz_pX(a, i, i+10, j);
  random_mat_zz_pX(b, i+10, i-10, j);

  t = GetTime();
  multiply_naive(c1, a, b);
  cout << GetTime()-t << " ";

  t = GetTime();
  multiply_evaluate(c2, a, b);
  cout << GetTime()-t << " ";

  if (c1 != c2)
    cout << "mismatch 1\n";

  zz_p::FFTInit(1);

  random_mat_zz_pX(a, i, i+10, j);
  random_mat_zz_pX(b, i+10, i-10, j);

  multiply_naive(c1, a, b);
  t = GetTime();
  multiply_evaluate(c2, a, b);
  cout << GetTime()-t << endl;

  if (c1 != c2)
    cout << "mismatch 2\n";
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
