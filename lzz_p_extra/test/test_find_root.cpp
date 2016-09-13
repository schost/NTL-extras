#include <NTL/lzz_p.h>

#include "lzz_p_extra.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/* if opt = 1, runs a check                                   */
/* else, runs timings                                         */
/*------------------------------------------------------------*/
void check(int opt){
   zz_p::init(9001);

   for (long i = 1; i < 1000; i += 1){

     zz_p a, b;

     if (opt == 1){
       find_root(a, i);
       find_root(b, i, a);
       for (long j = 1; j <= i; j++)
	 if (power(a, j) == to_zz_p(1) || power(b, j) == to_zz_p(1))
	   Error("wrong order in find_root");
     }
     else{
       cout << i << " ";

       double t;

       t = GetTime();
       for (long j = 0; j < 10000; j++){
	 find_root(a, i);
	 find_root(b, i, a);
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
