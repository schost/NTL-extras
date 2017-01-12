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

   for (long i = 10; i < 1100; i += 20){

     long NB = 1;
     if (i < 1000)
       NB = 5;
     if (i < 800)
       NB = 10;
     if (i < 700)
       NB = 20;
     if (i < 600)
       NB = 60;
     if (i < 500)
       NB = 125;
     if (i < 400)
       NB = 250;
     if (i < 300)
       NB = 500;
     if (i < 200)
       NB = 1000;
     if (i < 100)
       NB = 2000;
     if (i < 50)
       NB = 5000;

     Mat<zz_p> check, check_i;
     random(check, i, i);
     double t0 = GetTime();
     for (long j = 0; j < NB; j++)
       check_i = inv(check);
     t0 = GetTime() - t0;
     t0 = t0/NB;
     cout << i << " " << 0 << " " << t0 << endl;

     double t1 = 0;

     zz_p a = to_zz_p(9);
     for (long alpha = 1; alpha < i && t1 < 1.1*t0; alpha += 1){
       mat_zz_p A, B;
       random(A, i, alpha);
       random(B, i, alpha);
       lzz_p_cauchy_like_geometric M(A, B, to_zz_p(1), power(a, i), a);

       cout << i << " " << alpha << " ";
     
       t1 = GetTime();
       for (long k = 0; k < NB; k++){
	 lzz_p_cauchy_like_geometric Minv;
	 invert(Minv, M);
       }
       t1 = GetTime() - t1;
       t1 = t1/NB;
       cout << t1 << " ";
       
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
