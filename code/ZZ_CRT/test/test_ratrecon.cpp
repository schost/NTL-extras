#include <NTL/ZZ.h>
#include <NTL/vector.h>

#include "magma_output.h"
#include "ZZ_CRT.h"
#include "ratrecon.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/* does a multimod                                            */
/* check takes an extra argument, not used here               */
/*------------------------------------------------------------*/
void check(int opt){

  long sz = 10;
  long prime = 101;
  Vec<long> check;
  check.SetLength(3*sz);
  for (long i = 0; i < sz; i++)
    check[i] = InvMod(200 % prime, prime);
  for (long i = sz; i < 2*sz; i++)
    check[i] = 0;
  for (long i = 2*sz; i < 3*sz; i++)
    check[i] = InvMod(2000 % prime, prime);

  Vec<Vec<long>> ccheck;
  ccheck.SetLength(2);
  ccheck[0] = check;
  ccheck[1] = check;

  bool done = false;
  ZZ m = to_ZZ(100);
  while(done == false){
    m = NextPrime(m*10);
    Vec<ZZ> v;
    v.SetLength(3*sz);
    for (long i = 0; i < sz; i++)
      v[i] = InvMod(to_ZZ(200) % m, m);
    for (long i = sz; i < 2*sz; i++)
      v[i] = to_ZZ(0);
    for (long i = 2*sz; i < 3*sz; i++)
      v[i] = InvMod(to_ZZ(2000) % m, m);
    
    Vec<Vec<ZZ>> vv;
    vv.SetLength(2);
    vv[0] = v;
    vv[1] = v;

    Vec<ZZ> w;
    Vec<Vec<ZZ>> ww;
    ZZ den;
    // long s = ReconstructRational(ww, den, vv, m, ccheck, prime);
    long s = ReconstructRational(w, den, v, m, check, prime);
    if (s == 0)
      cout << "false with m=" << m << endl;
    else{
      cout << "true with m=" << m << endl;
      done = true;
      cout << den << endl;
      // cout << vv << endl;
      // cout << ww << endl;
      cout << w << endl;
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
