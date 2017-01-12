#ifndef ZZXY__H
#define ZZXY__H

#include <fstream>
#include <iostream>
#include <string.h>

#include <NTL/ZZX.h>
#include <NTL/vector.h>

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* a helper class for bivariate polynomials over ZZ           */
/* a ZZXY is simply a vector of ZZX                           */
/* with the convention f = sum_i coeffX[i](X) Y^i             */ 
/* minimal functionalities are provided                       */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

class ZZXY{
 public:
  Vec<ZZX> coeffX;

  /*------------------------------------------------------------*/
  /* total degree                                               */
  /*------------------------------------------------------------*/
  long tdeg() const;

  /*------------------------------------------------------------*/
  /* degree in Y                                                */
  /*------------------------------------------------------------*/
  long degY() const;

  /*------------------------------------------------------------*/
  /* degree in X                                                */
  /*------------------------------------------------------------*/
  long degX() const;

  /*------------------------------------------------------------*/
  /* naive evaluation algorithm to compute F(x,g(x)) mod x^t    */ 
  /*------------------------------------------------------------*/
  void eval(ZZX & val, const ZZX & g, long t);

  /*------------------------------------------------------------*/
  /* naive evaluation algorithm to compute F(x,g(x)/d) mod x^t  */ 
  /*------------------------------------------------------------*/
  void eval(ZZX & val, ZZ & dval, const ZZX & g, const ZZ & d, long t);

  /*------------------------------------------------------------*/
  /* finds g such that F(x,g(x)) = 0 mod x^t, g(0) = g0         */ 
  /*------------------------------------------------------------*/
  void series_solution(ZZX & g, ZZ& den_g, const ZZ & g0, long t);

  /*------------------------------------------------------------*/
  /* resizes the array of coefficients                          */
  /* to remove the trailing entries that are zero, if any       */
  /*------------------------------------------------------------*/
  void normalize();

  ZZXY(){}

  /*------------------------------------------------------------*/
  /* builds from a vector of ZZX                                */
  /*------------------------------------------------------------*/
  ZZXY(const Vec<ZZX> & coeff){
    coeffX = coeff;
    normalize();
  }

  /*------------------------------------------------------------*/
  /* builds from a file                                         */
  /* format: [[a0,a1,..,aN],...,[x0,x1,..,xM]]                  */
  /* where [a0,a1,..,aN] = coeff(f,Y^0) \in Z[X], etc.          */
  /*------------------------------------------------------------*/
  void read_from_file(const string& filename);

  ~ZZXY(){}
};

/*------------------------------------------------------------*/
/* I/O using NTL's representation (without commas!)           */
/*------------------------------------------------------------*/
istream& operator>>(istream& s, ZZXY& x);
ostream& operator<<(ostream& s, const ZZXY& a);

/*------------------------------------------------------------*/
/* initializes M<y,x>=QQ[y,x] with lex order y > x            */
/*------------------------------------------------------------*/
void magma_init_bi_QQ();

/*------------------------------------------------------------*/
/* prints a poly in M = QQ[y > x]                             */
/*------------------------------------------------------------*/
void magma_output(const ZZXY & a);

/*------------------------------------------------------------*/
/* assigns a poly in M = QQ[y > x]                            */
/*------------------------------------------------------------*/
void magma_assign(const ZZXY & a, const string & name);

/*------------------------------------------------------------*/
/* derivative in Y                                            */
/*------------------------------------------------------------*/
void diffY(ZZXY & dyF, const ZZXY & F);

/*------------------------------------------------------------*/
/* random poly with deg(F,x) <= dx, deg(F,y) <= dy            */
/* coefficients are in -2^{bound-1}..2^{bound-1}-1            */
/*------------------------------------------------------------*/
void random(ZZXY & F, long bound, long dx, long dy);

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* a class to represent the solutions of bivariate systems    */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
struct ZZ_bivariate_regular_chain{
  ZZX T1;
  ZZXY T2;
};

/*------------------------------------------------------------*/
/* prints a pair in M = QQ[y > x]                             */
/*------------------------------------------------------------*/
void magma_output(const ZZ_bivariate_regular_chain& T);

/*------------------------------------------------------------*/
/* assigns a pair in M = QQ[y > x]                            */
/*------------------------------------------------------------*/
void magma_assign(const ZZ_bivariate_regular_chain& T, const string& name);

/*------------------------------------------------------------*/
/* prints a vector of pairs in M = QQ[y > x]                  */
/*------------------------------------------------------------*/
void magma_output(const Vec<ZZ_bivariate_regular_chain>& T);

/*------------------------------------------------------------*/
/* assigns a vector of pairs in M = QQ[y > x]                 */
/*------------------------------------------------------------*/
void magma_assign(const Vec<ZZ_bivariate_regular_chain>& T, const string& name);

/*------------------------------------------------------------*/
/* computes a family of regular chains                        */
/* whose solutions is V(FZ, GZ)                               */
/* assumes V(FZ, GZ) finite, not empty                        */
/* otherwise, sols.length is set to zero                      */
/* the decomposition may not be the equiprojectable one       */
/*------------------------------------------------------------*/
void solve(Vec<ZZ_bivariate_regular_chain> & sols, const ZZXY & FZ, const ZZXY & GZ);

#endif
