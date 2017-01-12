#include <NTL/lzz_p.h>
#include <NTL/lzz_pX.h>
#include <assert.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
void check(int opt){

  zz_p::FFTInit(0);

  for (long i = 1; i < 1000; i+=1){
    Vec<zz_p> valf, valf_t;
    zz_pX f = random_zz_pX(i);
    zz_pX g = random_zz_pX(i);

    long NB = 40000;
    if (i < 100)
      NB = 100000;
    if (i < 10)
      NB = 500000;

    long n = 2*i;
    double t, u;
    zz_pX h;
    zz_pX_Multipoint_CTFT c;

    double t1, t2, t3, t4, t5;
    t1 = t2 = t3 = t4 = t5 = 0;

    u = GetTime();
    for (long k = 0; k < NB; k++)
      c = zz_pX_Multipoint_CTFT(n);
    t1 += GetTime()-u;

    t = GetTime();
    for (long k = 0; k < NB; k++){
      u = GetTime();
      c.evaluate(valf, f);
      t2 += GetTime()-u;

      u = GetTime();
      c.evaluate_t(valf_t, valf);
      t3 += GetTime()-u;

      u = GetTime();
      c.interpolate(f, valf);
      t4 += GetTime()-u;

      u = GetTime();
      c.interpolate_t(valf, valf_t);
      t5 += GetTime()-u;
    }
    t = (GetTime()-t)/NB;
    
    cout << i << " " << t1/NB << " " << t2/NB << " " << t3/NB << " " << t4/NB << " " << t5/NB;
    cout << endl;
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
