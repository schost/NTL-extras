#include "ZZ_p_block_sylvester.h"

NTL_CLIENT

/*----------------------------------------------------*/
/* only sets the type and output precision            */  
/*----------------------------------------------------*/
ZZ_p_block_sylvester::ZZ_p_block_sylvester(const Vec<long> &t, long p):
  type {t},
  prec {p}, 
  initialized {true} {
    max_of_type = 0;
    for (long i = 0; i < t.length(); i++)
      max_of_type = max(max_of_type, t[i]);
}

/*----------------------------------------------------*/
/* does nothing                                       */  
/*----------------------------------------------------*/
ZZ_p_block_sylvester::ZZ_p_block_sylvester(){
  initialized = false;
}

/*----------------------------------------------------*/
/* partitions the given vector into blocks            */
/* converts each block into polynomials               */
/*----------------------------------------------------*/
void ZZ_p_block_sylvester::create_lhs_list (Vec<ZZ_pX> & result, const Vec<ZZ_p> &v){
  const ZZ_p *vc = v.elts();
  result.SetLength(type.length());
  for (long i = 0; i < type.length(); i++){
    result[i].rep.SetLength(type[i] + 1);
    ZZ_p * cf = result[i].rep.elts();
    for (int j = 0; j < type[i] + 1; j++)
      cf[j] = vc[j];
    vc += type[i] + 1;
    result[i].normalize();
  }
}
