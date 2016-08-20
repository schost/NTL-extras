#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* does a multipoint evaluation                               */
/*------------------------------------------------------------*/
void check(int opt){

  long p = 1125899906842679;
  zz_p::init(p);

  if (opt == 1){
    for (long j = 1; j < 10; j++){
      Vec<zz_p> q;
      q.SetLength(j);
      for (long i = 0; i < j; i++){
	q[i] = random_zz_p();
      }
      
      {
	zz_pX_Multipoint * ev;
	zz_pX_Multipoint_General evQ(q);
	ev = &evQ;
	zz_pX f = random_zz_pX(2*j);
	Vec<zz_p> val;
	val.SetLength(j);
	ev->evaluate(val, f);
	for (long i = 0; i < j; i++)
	  cout << val[i]-eval(f, q[i]) << " ";
	cout << endl;
      }
      {
	zz_pX_Multipoint_General ev(q);
	zz_pX f = random_zz_pX(2*j);
	Vec<zz_p> val;
	val.SetLength(j);
	ev.evaluate(val, f);
	for (long i = 0; i < j; i++)
	  cout << val[i]-eval(f, q[i]) << " ";
	cout << endl;
      }
    }
  }
  else{
    double t;

    long j = 100;
    Vec<zz_p> q;
    q.SetLength(j);
    for (long i = 0; i < j; i++){
      q[i] = random_zz_p();
    }
    {
      zz_pX_Multipoint * ev;
      zz_pX_Multipoint_General evQ(q);
      ev = &evQ;
      zz_pX f = random_zz_pX(2*j);
      Vec<zz_p> val;
      val.SetLength(j);
      t = GetTime();
      ev->evaluate(val, f);
      cout << GetTime() - t << endl;
    }
    {
      zz_pX_Multipoint_General ev(q);
      zz_pX f = random_zz_pX(2*j);
      Vec<zz_p> val;
      val.SetLength(j);
      t = GetTime();
      ev.evaluate(val, f);
      cout << GetTime() - t << endl;
    }
  }
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
