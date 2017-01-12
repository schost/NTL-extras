#include <NTL/lzz_p.h>
#include <NTL/lzz_pX.h>
#include <assert.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
void check(int opt){

  zz_p::FFTInit(0);

  long start, end;
  if (opt == 1){
    start = 1;
    end = 700;
  }
  else{
    start = 1;
    end = 701;
  }
      

  for (long i = start; i < end; i+=1){
    Vec<zz_p> valf, valg, valh;
    zz_pX f = random_zz_pX(i);
    zz_pX g = random_zz_pX(i);

    long NB;
    if (opt == 1)
      NB = 1;
    else 
      NB = 4000;

    long n = 2*i;
    double t, t1;
    zz_pX h;
    zz_pX_Multipoint_CTFT c;


    t1 = GetTime();
    for (long k = 0; k < NB; k++)
      c = zz_pX_Multipoint_CTFT(n);
    t1 = GetTime()-t1;

    t = GetTime();
    for (long k = 0; k < NB; k++)
      c.mul(h, f, g);
    t = (GetTime()-t)/NB;
    
    if (opt == 1)
      cout << i << " " << (h-(f*g)) << endl;
    else{
      cout << i << " " << t << " ";
      t = GetTime();
      for (long k = 0; k < NB; k++)
	h = f*g;
      t = (GetTime()-t)/NB;
      cout << t << endl;
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
