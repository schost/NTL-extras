#include <vector>

#include <NTL/ZZ_pX.h>
#include <NTL/vector.h>

#include "ZZ_pX_extra.h"
#include "lzz_pX_CRT.h"
#include "ZZ_extra.h"
#include "magma_output.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/* prints a vector                                            */
/*------------------------------------------------------------*/
void magma_output(const Vec<ZZ_p> & v){
  if (v.length() == 0){
    cout << "[]";
    return;
  }
  cout << "[GF(" << ZZ_p::modulus() << ")|";
  for (long i = 0; i < v.length()-1; i++)
    cout << v[i] << ", ";
  cout << v[v.length()-1] << "]";
}

/*------------------------------------------------------------*/
/* assigns a vector to variable "name"                        */
/*------------------------------------------------------------*/
void magma_assign(const Vec<ZZ_p> & v, const string & name){
  cout << name << " := ";
  magma_output(v);
  cout << ";" << endl;
}

/*------------------------------------------------------------*/
/* prints a poly with indeterminate var                       */
/*------------------------------------------------------------*/
void magma_output(const ZZ_pX & v, const string & var){
  cout << "(Parent(" << var << ")!(0)";
  for (long i = 0; i <= deg(v); i++)
    cout << "+(" << coeff(v, i) << ")*" << var << "^" << i;
  cout << ")";
}

/*------------------------------------------------------------*/
/* prints a poly with indeterminate x, as a cast into U       */
/*------------------------------------------------------------*/
void magma_output(const ZZ_pX & v){
  cout << "U!(";
  magma_output(v.rep);
  cout << ")";
}
  
/*------------------------------------------------------------*/
/* assign a poly with indeterminate var to variable "name"    */
/*------------------------------------------------------------*/
void magma_assign(const ZZ_pX & v, const string & var, const string & name){
  cout << name << " := ";
  magma_output(v, var);
  cout << ";" << endl;
} 

/*------------------------------------------------------------*/
/* assign a poly with indeterminate x to variable "name"      */
/*------------------------------------------------------------*/
void magma_assign(const ZZ_pX & v, const string & name){
  magma_assign(v, "x", name);
}

/*------------------------------------------------------------*/
/* prints a vector of polys with indeterminate var            */
/*------------------------------------------------------------*/
void magma_output(const Vec<ZZ_pX> & v, const string & var){
  cout << "[";
  for (long i = 0; i < v.length(); i++){
    magma_output(v[i]);
    if (i < v.length()-1)
      cout << ", ";
  }
  cout << "]";
}

/*------------------------------------------------------------*/
/* prints a vector of polys with indeterminate x              */
/*------------------------------------------------------------*/
void magma_output(const Vec<ZZ_pX> & v){
  magma_output(v, "x");
}

/*------------------------------------------------------------*/
/* assign a vector of polys with indet var to variable "name" */
/*------------------------------------------------------------*/
void magma_assign(const Vec<ZZ_pX> & v, const string & var, const string & name){
  cout << name << ":=";
  magma_output(v, var);
  cout << ";\n";
}

/*------------------------------------------------------------*/
/* assign a vector of polys with indet x to variable "name"   */
/*------------------------------------------------------------*/
void magma_assign(const Vec<ZZ_pX> & v, const string & name){
  magma_assign(v, "x", name);
}

/*------------------------------------------------------------*/
/* maximum size of the entries of a                           */
/*------------------------------------------------------------*/
long size(const ZZ_pX& a){
  long d = 0;

  for (long j = 0; j <= deg(a); j++)
    d = max(d, size(coeff(a, j)._ZZ_p__rep));

  return d;
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* class ZZ_pX_poly_multiplier to multiply by a fixed argument*/
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* C = A*B                                                    */
/* assumes that deg(B) < n                                    */
/*------------------------------------------------------------*/
void ZZ_pX_poly_multiplier::mul(ZZ_pX& C, const ZZ_pX& B){
  zz_pPush push; 

  const ZZ_pFFTInfoT *FFTInfo = ZZ_p::GetFFTInfo();
  ZZ_pTmpSpaceT *TmpSpace = ZZ_p::GetTmpSpace();
  long nB = deg(B);

  Vec<zz_p> valC;
  valC.SetLength(n);

  const ZZ_p *xx = B.rep.elts();
  for (long j = 0; j <= nB; j++) {
    to_modular_rep(tmp, xx[j], FFTInfo, TmpSpace);
    for (long i = 0; i < nprimes; i++) 
      t[i][j] = tmp[i];
  }
  for (long j = nB+1; j < n; j++) 
    for (long i = 0; i < nprimes; i++) 
      t[i][j] = 0;
  
  for (long i = 0; i < nprimes; i++){
    ctxs[i].restore();
    Vec<zz_p> val;
    ctfts[i].evaluate(val, t[i].elts(), n);
    for (long j = 0; j < n; j++)
      valC[j] = valA[i][j] * val[j];
    zz_pX f;
    ctfts[i].interpolate(f, valC);
    long df = deg(f);
    zz_p * cf = f.rep.elts();
    for (long j = 0; j <= df; j++)
      t[i][j] = cf[j]._zz_p__rep;
    for (long j = df+1; j < n; j++)
      t[i][j]  = 0;
  }

  C.rep.SetLength(n);
  for (long j = 0; j < n; j++){
    for (long i = 0; i < nprimes; i++)
      tmp[i] = t[i][j];
    from_modular_rep(C.rep[j], tmp, FFTInfo, TmpSpace);
  }
  C.normalize();
}

/*------------------------------------------------------------*/
/* constructor, given an upper bound on the degree of args    */
/*------------------------------------------------------------*/
ZZ_pX_poly_multiplier::ZZ_pX_poly_multiplier(const ZZ_pX& A, long nB){

  zz_pPush push; 

  const ZZ_pFFTInfoT *FFTInfo = ZZ_p::GetFFTInfo();
  ZZ_pTmpSpaceT *TmpSpace = ZZ_p::GetTmpSpace();
  nprimes = FFTInfo->NumPrimes;
  
  ctxs.SetLength(nprimes);
  for (long i = 0; i < nprimes; i++)
    ctxs[i] = zz_pContext(INIT_FFT, i);

  long nA = deg(A);
  n = nA + nB;
  ctfts.SetLength(nprimes);
  for (long i = 0; i < nprimes; i++){
    ctxs[i].restore();
    ctfts[i] = zz_pX_Multipoint_CTFT(n);
  }

  valA.SetLength(nprimes);
  t.SetLength(nprimes);
  tmp.SetLength(nprimes);

  for (long i = 0; i < nprimes; i++){
    t[i].SetLength(n);
    valA[i].SetLength(n);
  }

  const ZZ_p *xx = A.rep.elts();
  for (long j = 0; j <= nA; j++) {
    to_modular_rep(tmp, xx[j], FFTInfo, TmpSpace);
    for (long i = 0; i < nprimes; i++) 
      t[i][j] = tmp[i];
  }
  for (long j = nA+1; j < n; j++) 
    for (long i = 0; i < nprimes; i++) 
      t[i][j] = 0;

  for (long i = 0; i < nprimes; i++){
    ctxs[i].restore();
    Vec<zz_p> val;
    ctfts[i].evaluate(val, t[i].elts(), n);
    for (long j = 0; j < n; j++)
      valA[i][j] = val[j];
  }
}
