#include <NTL/tools.h>
#include <time.h>
#include <cmath>
#include "bivariatemodularcomp.h"
#include "mat_ZZX.h"
#include "mat_ZZ_pX.h"
#include "mat_ZZ_pX_extra.h"
#include "ZZ_pX_extra.h"

using namespace std;
using namespace NTL;

void BivariateModularComp::init (const ZZ_pX &f, const Vec<long> &type_new, long prec_new){
  type = type_new;
  f_field = f;

  prec = prec_new;
  sqrtP = ceil(sqrt(type.length()));
  ZZ_pX running = ZZ_pX{1};
  Vec<ZZ_pX> fs;

  fs.SetLength(sqrtP);
  ZZ_pX_poly_multiplier multF(f, prec);
  for (long i = 0; i < sqrtP; i++){
    fs[i] = running;
    multF.mul(running, running);
    trunc(running, running, prec);
  }

  F_field = running;
  create_rhs_matrix(B, fs);
  initialized = true;
  mode = 0;
}

void BivariateModularComp::init(const Vec<ZZ_pX>& fs, const Vec<long> &type, long prec){
  sqrtP = ceil(sqrt(type.length()));
  if (fs.length() <= sqrtP) throw "not enough functions";
  this->type = type;
  this->prec = prec;
  fs_field = fs;
  initialized = true;
  mode = 1;
}

BivariateModularComp::BivariateModularComp(){}

BivariateModularComp::BivariateModularComp(const ZZ_pX& f, const Vec<long> &type, long prec) {
  init(f,type,prec);
}

BivariateModularComp::BivariateModularComp(const Vec<ZZ_pX>& fs, const Vec<long> &type, long prec){
  init(fs,type,prec);
}

void BivariateModularComp::create_lhs_list (Vec<ZZ_pX> & result, const Vec<ZZ_p> &v){
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

void BivariateModularComp::create_lhs_matrix(Mat<ZZ_pX> &A,const Vec<ZZ_pX> &v){
  A.SetDims(ceil((type.length()*1.0) / sqrtP), sqrtP);
  for (long i = 0; i < A.NumRows(); i++){
    for (long j = 0; j < A.NumCols(); j++){
      if(A.NumCols() * i + j < v.length()){
        A.put(i, j, v[A.NumCols() * i + j]);
      }
    }
  }
}

void BivariateModularComp::create_rhs_matrix(Mat<ZZ_pX> &B,const Vec<ZZ_pX> &v){
  B.SetDims(v.length(),type.length());
  for (long i = 0; i < B.NumRows(); i++){
    long acc = 0;
    for (long j = 0; j < B.NumCols(); j++){
      ZZ_pX p;
      for (int s = 0; s < type[j]+1; s++){
        SetCoeff(p,s,coeff(v[i],acc+s));
      }
      acc += type[j]+1;
      B.put(i,j,p);
    }
  }
}

void BivariateModularComp::deslice(Vec<ZZ_pX> &D, const Mat<ZZ_pX> &C){

  D.SetLength(C.NumRows());

  long s = 0;
  for (long i = 0; i < type.length(); i++)
    s += max(s, type[i]);

  for (long i = 0; i < C.NumRows(); i++){
    
    ZZ_pX v;
    v.rep.SetLength(prec + s + 1);
    ZZ_p* cv = v.rep.elts();
    const ZZ_pX * Ci = C[i].elts();
    const ZZ_p * cC = Ci[0].rep.elts();
    long old_len = Ci[0].rep.length();
    for (long j = 0; j < old_len; j++)
      cv[j] = cC[j];
    cv += type[0]+1;
    for (long a = 1; a < C.NumCols(); a++){
      const ZZ_p * cC = Ci[a].rep.elts();
      long new_len = Ci[a].rep.length();
      for (long j = 0; j < min(old_len-(type[a-1]+1), new_len); j++)
      	cv[j] += cC[j];
      if (old_len >= (type[a-1]+1))
	for (long j = old_len-(type[a-1]+1); j < new_len; j++)
	  cv[j] = cC[j];
      old_len = new_len;
      cv += type[a]+1;
    }

    v.normalize();
    D[i] = trunc(v, prec);
  }

  // D.SetLength(C.NumRows());
  // for (long i = 0; i < C.NumRows(); i++){
  //   D[i].rep.SetLength(prec);
  //   ZZ_p* cv = D[i].rep.elts();
  //   long acc = 0;
  //   for (long j = 0; j < C.NumCols(); j++){
  //     const ZZ_p* p = C[i][j].rep.elts();
  //     long dg = min(deg(C[i][j]) + 1, prec-acc);
  //     for (long s = 0; s < dg; s++)
  // 	cv[s] += p[s];
  //     cv += type[j]+1;
  //     acc += type[j]+1;
  //   }
  //   D[i].normalize();
  // }
}

Vec<ZZ_p> BivariateModularComp::mult_right_Horners(const Vec<ZZ_p> &rhs){
  if (!initialized) throw "must init first";
	if (mode != 0) throw "wrong mode";
  Vec<ZZ_pX> rhs_poly;
  create_lhs_list(rhs_poly, rhs);
  ZZ_pX result = rhs_poly[rhs_poly.length()-1];
  for (long i = rhs_poly.length()-2; i >= 0; i--){
    result = trunc((result * f_field), prec) + rhs_poly[i];
  }

  Vec<ZZ_p> v;
  v.SetLength(prec);
  for (long i = 0; i < prec; i++)
    v[i] = coeff(result, i);
  return v;
}

Vec<ZZ_p> BivariateModularComp::mult_right_naive(const Vec<ZZ_p> &rhs){
	if (!initialized) throw "must init first";
	Vec<ZZ_pX> rhs_poly;
	create_lhs_list(rhs_poly,rhs);
	ZZ_pX result;
	for (long i = 0; i < fs_field.length(); i++)
		result += fs_field[i] * rhs_poly[i];
	Vec<ZZ_p> v;
	v.SetLength(prec);
	for (long i = 0; i < prec; i++)
		v[i] = coeff(result, i);
	return v;
}

Vec<ZZ_p> BivariateModularComp::mult_right(const Vec<ZZ_p> &rhs){
	if (mode == 0) return mult_right_comp(rhs);
	return mult_right_naive(rhs);
}

Vec<ZZ_p> BivariateModularComp::mult_right_comp(const Vec<ZZ_p> &rhs){
  if (!initialized) throw "must init first";
  if (mode != 0) throw "wrong mode";
  Mat<ZZ_pX> A0, A;
  Vec<ZZ_pX> rhs_poly;

  create_lhs_list(rhs_poly, rhs);
  create_lhs_matrix(A0, rhs_poly);
  mul_CRT_CTFT(A, A0, B);
  Vec<ZZ_pX> B1; // TODO: change this name!

  deslice(B1, A);

  ZZ_pX p;
  p = B1[B1.length() - 1];
  ZZ_pX_poly_multiplier multF(F_field, prec);
  for (long i = B1.length()-2; i >= 0; i--){
    multF.mul(p, p);
    p = trunc(p, prec) + B1[i];
  }

  Vec<ZZ_p> v;
  v.SetLength(prec);
  for (long i = 0; i < prec; i++)
    v[i] = coeff(p, i);

  return v;
}




