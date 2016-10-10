#include <NTL/lzz_p.h>
#include <NTL/lzz_pX.h>
#include <assert.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
void check(int opt){

  zz_p::FFTInit(0);

  for (long i = 2; i < 1000; i+=1){
    zz_pX f;
    SetCoeff(f, 1, 1);

    zz_pX_Multipoint_CTFT c(i);
    
    Vec<zz_p> val;
    c.evaluate(val, f);
    cout << i << endl;
    assert (val == c.pts);
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
