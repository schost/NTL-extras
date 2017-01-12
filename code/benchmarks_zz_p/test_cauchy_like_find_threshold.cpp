#include <NTL/vec_lzz_p.h>
#include <assert.h>

#include "vec_lzz_p_extra.h"
#include "lzz_p_cauchy_geometric.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* if opt = 1, runs FFT0                                      */
/* else, runs mod 65537                                       */
/*------------------------------------------------------------*/
void check(int opt){
  if (opt == 1)
    zz_p::FFTInit(0);
  else
    zz_p::UserFFTInit(65537);
  
  for (long alpha = 2; alpha < 60; alpha += 4){
    for (long i = 1000; i < 4000; i += 1000){    
      zz_p a = to_zz_p(9);
      long j = i;
      
      mat_zz_p A, B;
      random(A, i, alpha);
      random(B, j, alpha);
      lzz_p_cauchy_like_geometric M(A, B, to_zz_p(1), power(a, i), a);
       
      cout << alpha << " " << i << " ";
      for (long thresh = 100; thresh < min(i, 1000); thresh += 30){
	double t;
	t = GetTime();
	lzz_p_cauchy_like_geometric Minv;
	invert_fast(Minv, M, thresh);
	t = GetTime() - t;
	cout << thresh << " " << t << " ";
      }
      cout << endl;
    }
  }
}

/*------------------------------------------------------------*/
/* main just calls check()                                    */
/* if not argument is given, runs timings                     */
/* if the argument 1 is given, runs check                     */
/*------------------------------------------------------------*/
int main(int argc, char** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);

  return 0;
}
