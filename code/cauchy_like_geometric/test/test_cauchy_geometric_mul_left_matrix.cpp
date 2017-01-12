#include <NTL/vec_lzz_p.h>
#include <assert.h>

#include "vec_lzz_p_extra.h"
#include "lzz_p_cauchy_geometric.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* if opt = 1, runs a check                                   */
/* else, runs timings                                         */
/*------------------------------------------------------------*/
void check(int opt){
   zz_p::FFTInit(0);

   for (long i = 1; i < 10000; i += 1){

     zz_p a = random_zz_p();

     if (opt == 1){
       lzz_p_cauchy_geometric C;
       Mat<zz_p> M;
       Mat<zz_p> input, output, check;

       long j = i+19;
       long k = i+9;

       random(input, i, k);
       C = lzz_p_cauchy_geometric(to_zz_p(1), power(a, i), a, i, j);
       C.to_dense(M);
       C.mul_left(output, input);
       check = transpose(M)*input;
       assert (check == output);

       random(input, j, k);
       C = lzz_p_cauchy_geometric(to_zz_p(1), power(a, i+19), a, j, i);
       C.to_dense(M);
       C.mul_left(output, input);
       check = transpose(M)*input;
       assert (check == output);

       cout << i << endl;
     }
     else{
       cout << i << " ";

       double t;

       t = GetTime();
       for (long j = 0; j < 10000; j++){
	 lzz_p_cauchy_geometric C(to_zz_p(1), power(a, i), a, i, i+19);
	 Mat<zz_p> M;
	 C.to_dense(M);
       }
       t = GetTime() - t;
       cout << t << " ";

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
