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

   for (long i = 10; i < 1000; i += 20){
     for (long alpha = 4; alpha < min(i, 60); alpha += 1){
       zz_p a = to_zz_p(9);
       long j = i;
       mat_zz_p A, B;
       random(A, i, alpha);
       random(B, j, alpha);
       lzz_p_cauchy_like_geometric M(A, B, to_zz_p(1), power(a, i), a);

       cout << i << " " << alpha << " ";
       double t;
       long NB = 1;
       if (i < 200)
	 NB = 200;
       if (i < 100)
	 NB = 500;
       if (i < 50)
	 NB = 10000;

       t = GetTime();
       for (long k = 0; k < NB; k++){
	 lzz_p_cauchy_like_geometric Minv;
	 invert_block(Minv, M);
       }
       t = GetTime() - t;
       cout << t/NB << " ";
        
       t = GetTime();
       for (long k = 0; k < NB; k++){
	 lzz_p_cauchy_like_geometric Minv;
	 invert_direct(Minv, M);
       }
       t = GetTime() - t;
       cout << t/NB << " ";

       // t = GetTime();
       // for (long k = 0; k < 1; k++){
       // 	 Mat<zz_p> Zinv;
       // 	 inv(Zinv, Z);
       // }
       // t = GetTime() - t;
       // cout << t << " ";

       // vec_zz_p in, out;
       // zz_p d;
       // random(in, i);
       // t = GetTime();
       // for (long k = 0; k < 1; k++){
       // 	 solve(d, Z, out, in);
       // }
       // t = GetTime() - t;
       // cout << t << " ";
       
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
