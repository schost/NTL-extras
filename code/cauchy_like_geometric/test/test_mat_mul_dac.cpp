#include <NTL/vec_lzz_p.h>
#include <assert.h>

#include "vec_lzz_p_extra.h"
#include "lzz_pX_cauchy_geometric_special.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* if opt = 1, runs a check                                   */
/* else, runs timings                                         */
/*------------------------------------------------------------*/
void check(int opt){
   zz_p::FFTInit(0);
   zz_p::UserFFTInit(65537);

   for (long n = 128; n < 16384; n = 2*n){
     for (long alpha = 2; alpha < n/4; alpha = 2*alpha){
       double t;
       Vec<Vec<zz_p>> A, Acheck, C, G, H;
       A.SetLength(n);
       C.SetLength(n);
       G.SetLength(n);
       H.SetLength(n);
       
       for (long i = 0; i < n; i++){
	 A[i].SetLength(alpha);
	 C[i].SetLength(alpha);
	 G[i].SetLength(alpha);
	 H[i].SetLength(alpha);
	 
	 for (long j = 0; j < alpha; j++){
	   A[i][j] = random_zz_p();
	   C[i][j] = random_zz_p();
	   G[i][j] = random_zz_p();
	   H[i][j] = random_zz_p();
	 }
       }
       
       cout << alpha << " " << n << " ";
       t = GetTime();
       mat_mul_sigmaLU(A, G, H, C);
       cout << GetTime()-t << " ";

       t = GetTime();
       mat_mul_sigmaLU_FFT(Acheck, G, H, C);
       cout << GetTime()-t << " ";
       cout << (Acheck == A) << " ";

       // zz_p a = random_zz_p();
       // cauchy_geometric_special Cc(to_zz_p(1), power(a, n), a, n, n);
       // mat_zz_p Aa, Bb;
       // random(Aa, n, alpha);
       // random(Bb, n, alpha);
       // cauchy_like_geometric_special M(Aa, Bb, Cc);
       // mat_zz_p Z;

       // vec_zz_p in, out;
       // random(in, n);
     
       // t = GetTime();
       // for (long j = 0; j < alpha; j++)
       // 	 mul_right(out, M, in);
       // t = GetTime() - t;
       // cout << t << " ";

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
