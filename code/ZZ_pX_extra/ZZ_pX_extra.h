#ifndef __ZZ_PX_EXTRA__H
#define __ZZ_PX_EXTRA__H

#include <NTL/ZZ_pX.h>
#include <NTL/lzz_pX.h>
#include <NTL/vector.h>


#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* prints a vector                                            */
/*------------------------------------------------------------*/
void magma_output(const Vec<ZZ_p> & v);

/*------------------------------------------------------------*/
/* assigns a vector to variable "name"                        */
/*------------------------------------------------------------*/
void magma_assign(const Vec<ZZ_p> & v, const string & name);

/*------------------------------------------------------------*/
/* prints a poly with indeterminate var                       */
/*------------------------------------------------------------*/
void magma_output(const ZZ_pX & v, const string & var);

/*------------------------------------------------------------*/
/* prints a poly with indeterminate x, as a cast into U       */
/*------------------------------------------------------------*/
void magma_output(const ZZ_pX & v);

/*------------------------------------------------------------*/
/* assign a poly with indeterminate var to variable "name"    */
/*------------------------------------------------------------*/
void magma_assign(const ZZ_pX & v, const string & var, const string & name);

/*------------------------------------------------------------*/
/* assign a poly with indeterminate x to variable "name"      */
/*------------------------------------------------------------*/
void magma_assign(const ZZ_pX & v, const string & name);

/*------------------------------------------------------------*/
/* prints a vector of polys with indeterminate var            */
/*------------------------------------------------------------*/
void magma_output(const Vec<ZZ_pX> & v, const string & var);

/*------------------------------------------------------------*/
/* prints a vector of polys with indeterminate x              */
/*------------------------------------------------------------*/
void magma_output(const Vec<ZZ_pX> & v);

/*------------------------------------------------------------*/
/* assign a vector of polys with indet var to variable "name" */
/*------------------------------------------------------------*/
void magma_assign(const Vec<ZZ_pX> & v, const string & var, const string & name);

/*------------------------------------------------------------*/
/* assign a vector of polys with indet x to variable "name"   */
/*------------------------------------------------------------*/
void magma_assign(const Vec<ZZ_pX> & v, const string & name);

/*------------------------------------------------------------*/
/* maximum size of the entries of a                           */
/*------------------------------------------------------------*/
long size(const ZZ_pX& a);

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* class to multiply by a fixed argument                      */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

class ZZ_pX_poly_multiplier{
 public:
  /*------------------------------------------------------------*/
  /* C = A*B                                                    */
  /* assumes that deg(B) < n                                    */
  /*------------------------------------------------------------*/
  void mul(ZZ_pX& C, const ZZ_pX& B);

  /*------------------------------------------------------------*/
  /* constructor, given an upper bound on the degree of args    */
  /*------------------------------------------------------------*/
  ZZ_pX_poly_multiplier(const ZZ_pX& A, long nB);

 private:
  Vec<zz_pX_Multipoint_CTFT> ctfts;
  Vec<zz_pContext> ctxs;
  Vec<Vec<zz_p>> valA;
  Vec<Vec<long>> t;
  Vec<long> tmp;

  long nprimes;
  long n;  // maximum possible length of the product
};


#endif
