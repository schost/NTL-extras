#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* does a multipoint evaluation                               */
/*------------------------------------------------------------*/
void check(int opt){

  zz_p::FFTInit(0);

  
  for (long j = 1; j < 10; j++){
    Vec<zz_p> q;
    q.SetLength(j);
    for (long i = 0; i < j; i++)
      q[i] = random_zz_p();
    
    zz_pX_Multipoint_General ev(q);
    for (long i = 0; i < j; i++){
      zz_p pt;
      ev.point(pt, i);
      cout << pt-q[i] << " ";
    }
    cout << endl;
  }

  for (long j = 1; j < 10; j++){
    zz_p a = random_zz_p();
    zz_pX_Multipoint_Geometric ev(a, j, 22);
    for (long i = 0; i < j; i++){
      zz_p pt;
      ev.point(pt, i);
      cout << pt-power(a, 2*(i+22)) << " ";
    }
    cout << endl;
  }

  for (long j = 1; j < 10; j++){
    zz_p a = random_zz_p();
    zz_p b = random_zz_p();
    zz_pX_Multipoint_Geometric ev(a, j, b);
    for (long i = 0; i < j; i++){
      zz_p pt;
      ev.point(pt, i);
      cout << pt - b*power(a, 2*i) << " ";
    }
    cout << endl;
  }

  for (long j = 1; j < 20; j=2*j){
    zz_pX_Multipoint_FFT ev(j);
    Vec<zz_p> pts;
    zz_pX X;
    X = 0;
    SetCoeff(X, 1, 1);
    ev.evaluate(pts, X);
    for (long i = 0; i < j; i++){
      zz_p pt;
      ev.point(pt, i);
      cout << pt - pts[i] << " ";
    }
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
