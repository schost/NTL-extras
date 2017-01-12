#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_CRT.h"
#include "magma_output.h"
#include "lzz_pXY.h"
#include "lzz_pEX_augmented.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* check takes an extra argument, not used here               */
/*------------------------------------------------------------*/
void check(int opt){

  long p = 939000522306748417;
  zz_p::init(p);

  zz_pX c0,c1,c2;
  c0.rep.append(to_zz_p(145995537835310784));
  c0.rep.append(to_zz_p(760812726164890992));
  c0.rep.append(to_zz_p(179994688516630272));
  c0.rep.append(to_zz_p(826858574696249967));

  c1.rep.append(to_zz_p(754237903198973686));
  c1.rep.append(to_zz_p(643208041538984805));
  c1.rep.append(to_zz_p(469280123987029674));
  c1.rep.append(to_zz_p(722575577799590066));

  c2.rep.append(to_zz_p(754237903198973686));
  c2.rep.append(to_zz_p(643208041538984805));
  c2.rep.append(to_zz_p(469280123987029674));

  Vec<zz_pX> C;
  C.append(c0);
  C.append(c1);
  C.append(c2);

  zz_pX d0, d1, d2;
  d0.rep.append(to_zz_p(884603674384760447));
  d0.rep.append(to_zz_p(864814845505761066));
  d0.rep.append(to_zz_p(355195820442067657)); 
  d0.rep.append(to_zz_p(68719164061005511));

  d1.rep.append(to_zz_p(435732182695290395));
  d1.rep.append(to_zz_p(368654333399044787));
  d1.rep.append(to_zz_p(118853739087555615));
  d1.rep.append(to_zz_p(38924617328406598));

  d2.rep.append(to_zz_p(435732182695290395));
  d2.rep.append(to_zz_p(368654333399044787));
  d2.rep.append(to_zz_p(118853739087555615));
  
  Vec<zz_pX> D;
  D.append(d0);
  D.append(d1);
  D.append(d2);

  zz_pXY pC(C), pD(D), G;
   
  GCD(G, pC, pD);
  magma_init();
  magma_init_bi();
  magma_assign_bi(G, "G");
}  

int main(int argc, char ** argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
