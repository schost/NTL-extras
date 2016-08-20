#ifndef LZZ_PXY__H
#define LZZ_PXY__H

#include <NTL/ZZX.h>
#include <NTL/lzz_pX.h>
#include <NTL/lzz_pEX.h>
#include <NTL/vector.h>

#include "lzz_pX_CRT.h"
#include "ZZXY.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* a helper class for bivariate polynomials over zz_p         */
/* a zz_pXY is simply a vector of zz_pX                       */
/* with the convention f = sum_i coeffX[i](X) Y^i             */ 
/* minimal functionalities are provided                       */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

class zz_pXY{
 public:
  Vec<zz_pX> coeffX;

  /*------------------------------------------------------------*/
  /* total degree                                               */
  /*------------------------------------------------------------*/
  long tdeg() const;

  /*------------------------------------------------------------*/
  /* degree in X                                                */
  /*------------------------------------------------------------*/
  long degX() const;

  /*------------------------------------------------------------*/
  /* degree in Y                                                */
  /*------------------------------------------------------------*/
  long degY() const;

  /*------------------------------------------------------------*/
  /* resizes the array of coefficients                          */
  /* to remove the trailing entries that are zero, if any       */
  /*------------------------------------------------------------*/
  void normalize();

  zz_pXY(){}

  /*------------------------------------------------------------*/
  /* builds from a vector of zz_pX                              */
  /*------------------------------------------------------------*/
  zz_pXY(const Vec<zz_pX> & coeff){
    coeffX = coeff;
    normalize();
  }

  ~zz_pXY(){}
};

/*------------------------------------------------------------*/
/* I/O using NTL's representation (without commas!)           */
/*------------------------------------------------------------*/
ostream& operator<<(ostream& s, const zz_pXY& a);

/*------------------------------------------------------------*/
/* initializes V<y>=GF(p)[x][y]                               */
/*------------------------------------------------------------*/
void magma_init_Y();

/*------------------------------------------------------------*/
/* initializes M<y,x>=GF(p)[y,x] with lex order y > x         */
/*------------------------------------------------------------*/
void magma_init_bi();

/*------------------------------------------------------------*/
/* prints a poly in V = GF(p)[x][y]                           */
/*------------------------------------------------------------*/
void magma_output(const zz_pXY & a);

/*------------------------------------------------------------*/
/* assigns a poly in V = GF(p)[x][y]                          */
/*------------------------------------------------------------*/
void magma_assign(const zz_pXY & a, const string & n);

/*------------------------------------------------------------*/
/* prints a poly in M = GF(p)[y > x]                          */
/*------------------------------------------------------------*/
void magma_output_bi(const zz_pXY & a);

/*------------------------------------------------------------*/
/* assigns a poly in M = GF(p)[y > x]                         */
/*------------------------------------------------------------*/
void magma_assign_bi(const zz_pXY & a, const string & name);

/*------------------------------------------------------------*/
/* random element f with deg(f,x) < dx and deg(f,y) < dy      */
/*------------------------------------------------------------*/
void random_zz_pXY(zz_pXY & f, long dx, long dy);

/*------------------------------------------------------------*/
/* multipoint evaluation with respect to X                    */
/* result is a vector of polynomials                          */
/*------------------------------------------------------------*/
void evaluate(Vec<zz_pX> & values, const zz_pXY & f, const zz_pX_Multipoint & ev);

/*------------------------------------------------------------*/
/* computes the resultant of f and g                          */
/* together with the subresultant of degree 1, if possible    */
/*------------------------------------------------------------*/
void resultant(zz_pX & res, zz_pX& subres0, zz_pX& subres1, const zz_pXY& f, const zz_pXY& g);
/* void resultant(zz_pX & res, const zz_pXY& f, const zz_pXY& g); */

/*------------------------------------------------------------*/
/* derivative in Y                                            */
/*------------------------------------------------------------*/
void diffY(zz_pXY & dyF, const zz_pXY & F);

/*------------------------------------------------------------*/
/* conversions                                                */
/*------------------------------------------------------------*/
void conv(zz_pEX & fEX, const zz_pXY & f);
void conv(zz_pXY & f, const zz_pEX & fEX);
void conv(zz_pXY & f, const ZZXY & fZ);

/*------------------------------------------------------------*/
/* reduces F(x,y) modulo (R(x), y-S(x))                       */
/*------------------------------------------------------------*/
void reduce_naive(zz_pX & rem, const zz_pX& R, const zz_pX& S, const zz_pXY & F);
void reduce(zz_pX & rem, const zz_pX& R, const zz_pX& S, const zz_pXY & F);


#endif

