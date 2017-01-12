#include <NTL/lzz_pX.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <assert.h>

#include "vec_lzz_p_extra.h"
#include "lzz_p_toeplitz_like.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* initializes a toeplitz_like; test reconstruction           */
/* check takes an extra argument, not used here               */
/*------------------------------------------------------------*/
void check(int opt){

  zz_p::UserFFTInit(65537);

  for (long alpha = 10; alpha < 16; alpha ++)
    for (long m = 10; m < 30; m++)
      for (long n = m-3; n < m+9; n++){
	Mat<zz_p> G, H;
	random(G, m, alpha);
	random(H, n, alpha);
	lzz_p_toeplitz_like T(G, H);
	Mat<zz_p> M;
	T.to_dense(M);

	Mat<zz_p> input, output, check;
	random(input, n, alpha+2);
	mul(check, M, input);
	T.mul_right_dac(output, input);
	assert (output == check);
	cout << m << " " << n << " " << alpha << endl;
	random(input, n, alpha-2);
	mul(check, M, input);
	T.mul_right_dac(output, input);
	assert (output == check);
	cout << m << " " << n << " " << alpha << endl;
      }
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
