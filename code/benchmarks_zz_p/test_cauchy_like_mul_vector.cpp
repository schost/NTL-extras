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
  for (long i = 10; i < 4000; i += 20){

     long NB = 100;
     if (i < 3000)
       NB = 300;
     if (i < 2000)
       NB = 1000;
     if (i < 1000)
       NB = 5000;
     if (i < 500)
       NB = 10000;
     if (i < 100)
       NB = 100000;

     Mat<zz_p> check;
     Vec<zz_p> in, out;
     random(check, i, i);
     random(in, i);
     double t0 = GetTime();
     for (long j = 0; j < NB; j++)
       out = check*in;
     t0 = GetTime() - t0;
     t0 = t0/NB;
     cout << i << " " << 0 << " " << t0 << endl;

     double t = 0;
     zz_p a = to_zz_p(9);
     for (long alpha = 1; alpha < i && t < 1.1*t0; alpha += 1){
       mat_zz_p A, B;
       random(A, i, alpha);
       random(B, i, alpha);
       lzz_p_cauchy_like_geometric M(A, B, to_zz_p(1), power(a, i), a);

       Vec<zz_p> in, out;
       random(in, i);
       
       cout << i << " " << alpha << " ";

       t = GetTime();
       for (long j = 0; j < NB; j++)
	 M.mul_right(out, in);
       t = GetTime() - t;
       t = t/NB;
       cout << t << " ";
              
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
