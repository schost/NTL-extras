#include <NTL/lzz_pX.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <assert.h>

#include "vec_lzz_p_extra.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* initializes a zz_pX_Multipoint                             */
/* check takes an extra argument, not used here               */
/*------------------------------------------------------------*/
void check(int opt){

  long p = 1125899906842679;
  zz_p::init(p);

  for (long j = 10; j < 1100; j++){
    cout << j << " ";
    zz_p a = random_zz_p();
    zz_p c = random_zz_p();

    Mat<zz_p> M;
    M.SetDims(j, j);
    for (long i = 0; i < j; i++)
      for (long k = 0; k < j; k++)
	M[i][k] = power(c*power(a, 2*i), k);

    Vec<zz_p> input, output, check;
    random(input, j);

    zz_pX_Multipoint_Geometric evG(a, j, c);
    mul(output, input, M);
    evG.inverse_mul_left(check, output);

    assert (check == input);
    cout << endl;
  }
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
