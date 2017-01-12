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

  for (long j = 100; j < 110; j++){
    long nb = j*j;
    Vec<zz_p> q;
    zz_p a = random_zz_p();
    long i0 = 246;
    q.SetLength(nb);
    for (long i = 0; i < nb; i++){
      q[i] = power(a,2*(i+i0));
    }

    zz_pX f = random_zz_pX(j);
    Vec<zz_p> valQ, valG;

    double t;

    zz_pX_Multipoint * ev;
    zz_pX_Multipoint_General evQ(q);
    ev = &evQ;

    cout << j << " ";

    t = GetTime();
    for (long i = 0; i < 10; i++)
      ev->evaluate(valQ, f);
    cout << GetTime()-t << " ";

    zz_pX_Multipoint_Geometric evG(a, j, nb, i0);
    ev = &evG;
    t = GetTime();
    for (long i = 0; i < 10; i++)
      ev->evaluate(valG, f);
    cout << GetTime()-t << " ";
    if (valG != valQ) 
      cout << j << endl;

    Mat<zz_p> mata, matb, matc;
    mata.SetDims(j/2, j/2);
    matb.SetDims(j/2, j);
    for (long i = 0; i < j/2; i++)
      for (long k = 0; k < j/2; k++)
	mata[i][k] = random_zz_p();
    for (long i = 0; i < j/2; i++)
      for (long k = 0; k < j; k++)
	matb[i][k] = random_zz_p();

    t = GetTime();
    for (long i = 0; i < 10; i++){
      for (long k = 0; k < j; k++){
	mul(matc, mata, matb);
	mul(matc, mata, matb);
      }
    }
    t = (GetTime()-t)/j;
    cout << t << " ";


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
