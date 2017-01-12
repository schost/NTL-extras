#include <NTL/lzz_pX.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* initializes a zz_pX_Multipoint                             */
/* check takes an extra argument, not used here               */
/*------------------------------------------------------------*/
void check(int opt){

  long p = 1125899906842679;
  zz_p::init(p);

  for (long j = 30; j < 100; j++){
    long nb = j*j;
    Vec<zz_p> q;
    zz_p a = random_zz_p();
    long i0 = 246;
    q.SetLength(nb);
    for (long i = 0; i < nb; i++){
      q[i] = power(a,2*(i+i0));
    }

    Vec<zz_pX> f;
    f.SetLength(j);
    for (long k = 0; k < j; k++)
      f[k] = random_zz_pX(j);
    Vec<Vec<zz_p>> valQ, valG;

    double t;

    zz_pX_Multipoint * ev;
    zz_pX_Multipoint_General evQ(q);
    ev = &evQ;

    cout << j << " ";

    t = GetTime();
    for (long i = 0; i < 1; i++)
      ev->evaluate(valQ, f);
    cout << GetTime()-t << " ";

    t = GetTime();
    zz_pX_Multipoint_Geometric evG(a, j, nb, i0);
    cout << GetTime()-t << " ";
    ev = &evG;

    t = GetTime();
    for (long i = 0; i < 1; i++)
      ev->evaluate(valG, f);
    cout << GetTime()-t << " ";
    if (valG != valQ) 
      cout << j << endl;
    
    cout << endl;


  }
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
