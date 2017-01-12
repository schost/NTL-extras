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

   for (long i = 100; i < 5000; i += 50){
     for (long alpha = 2; alpha < 60; alpha += 2){
       zz_p a = to_zz_p(9);
       long j = i+10;

       mat_zz_p A, B;
       random(A, i, alpha);
       random(B, j, alpha);
       lzz_p_cauchy_like_geometric M(A, B, to_zz_p(1), power(a, i), a);

       Mat<zz_p> in, out;
       random(in, j, alpha);
       
       cout << i << " " << alpha << " ";
       double t;
       
       long NB = 1;
       t = GetTime();
       for (long j = 0; j < NB; j++)
	 M.mul_right(out, in);
       t = GetTime() - t;
       cout << t/NB << " ";

       t = GetTime();
       for (long j = 0; j < NB; j++)
	 M.mul_right_sigma_UL(out, in);
       t = GetTime() - t;
       cout << t/NB << " ";
              
       cout << endl;
     }
     cout << endl;
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
