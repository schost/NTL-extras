#include <NTL/lzz_p.h>
#include <NTL/lzz_pX.h>
#include <assert.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
void check(int opt){

  zz_p::FFTInit(0);

  while(1){
    for (long i = 1; i < 1000; i+=1){
      Vec<zz_p> val;
      zz_pX f = random_zz_pX(i);
      zz_pX_Multipoint_CTFT c(i);
      c.evaluate(val, f);
      
      zz_pX g;
      c.interpolate(g, val);
      cout << i << " " << f-g << endl;
    }
  }
}

/*------------------------------------------------------------*/
/* main just calls check()                                    */
/* if not argument is given, runs timings                     */
/* if the argument 1 is given, runs check                     */
/*------------------------------------------------------------*/
int main(int argc, char **argv){

  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);

  return 0;
}
