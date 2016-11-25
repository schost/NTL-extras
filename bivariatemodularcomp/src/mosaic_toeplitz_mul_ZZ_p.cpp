#include "mosaic_toeplitz_mul_ZZ_p.h"
using namespace NTL;

mosaic_toeplitz_mul_ZZ_p::mosaic_toeplitz_mul_ZZ_p(const Vec<long> &t,long p):type{t},prec{p}, initialized{true}{}

mosaic_toeplitz_mul_ZZ_p::mosaic_toeplitz_mul_ZZ_p(){initialized=false;}

void mosaic_toeplitz_mul_ZZ_p::create_lhs_list (Vec<ZZ_pX> & result, const Vec<ZZ_p> &v){
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
