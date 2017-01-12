#ifndef LZZ_PX_AUGMENTED__H
#define LZZ_PX_AUGMENTED__H

#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
#include <algorithm>

#include <NTL/ZZX.h>
#include <NTL/lzz_pX.h>
#include <NTL/lzz_pEX.h>
#include <NTL/vector.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/* helper function for splitting modulus into two factors:    */
/* one "mod_unit" where elt is a unit (we compute its inverse)*/
/* one "mod_zero" where elt is zero                           */
/* assumes modulus to be squarefree                           */
/*------------------------------------------------------------*/
void split(zz_pX & mod_zero, zz_pX & mod_unit, zz_pX & inverse, const zz_pX & elt, const zz_pX & modulus);

/*------------------------------------------------------------*/
/* Let the "degree" of a vector of zz_p be the index of the   */
/* last non-zero entry.                                       */
/* This function splits "modulus" into factors f[i]           */
/* such that for any two roots x,x' of f[i],                  */
/*  deg(elts(x))=deg(elts(x'))=:d[i]                           */
/* and d[i] != d[i'] for i != i'                              */
/* assumes modulus to be squarefree                           */
/*------------------------------------------------------------*/
void split_for_regularize(Vec<zz_pX> & factor_modulus, const Vec<zz_pX> & elts, const zz_pX & modulus);


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* a class similar to zz_pEX, with a copy of a zz_pEContext   */
/* minimal functionalities are provided                       */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class zz_pEX_augmented{
 public:
  zz_pX T1X; // keep a copy; cannot access it from T1 easily without switching moduli?
  zz_pEContext T1;
  zz_pEX T2;

  zz_pEX_augmented(){}

  zz_pEX_augmented(const zz_pEX & t2){
    T1X = zz_pE::modulus();
    T1 = zz_pEContext(T1X);
    T2 = t2;
  }

  zz_pEX_augmented(const zz_pX & t1, const zz_pXY & t2){
    T1X = t1;
    T1 = zz_pEContext(t1);
    zz_pEPush push(T1); 
    conv(T2, t2);
  }

  zz_pEX_augmented(const zz_pX & t1, const Vec<zz_pX> & t2){
    T1X = t1;
    T1 = zz_pEContext(t1);
    zz_pEPush push(T1); 

    Vec<zz_pE> t2E;
    conv(t2E, t2);
    conv(T2, t2E);
  }

  ~zz_pEX_augmented(){}
};

/*------------------------------------------------------------*/
/* outputs as a sequence of M<y,x>                            */
/*------------------------------------------------------------*/
void magma_outputs(const zz_pEX_augmented & F);

/*------------------------------------------------------------*/
/* assigns to "name" as a sequence of M<y,x>                  */
/*------------------------------------------------------------*/
void magma_assign(const zz_pEX_augmented & F, const string & name);

/*------------------------------------------------------------*/
/* random element with deg(T1,x) = dX and deg(T2,y) = dY      */
/*------------------------------------------------------------*/
void random_monic_zz_pEX_augmented(zz_pEX_augmented & res, long dX, long dY);

/*------------------------------------------------------------*/
/* replaces T1 by new_mod and reduces T2 by it                */
/*------------------------------------------------------------*/
void mod(zz_pEX_augmented & reduced, const zz_pEX_augmented & F, const zz_pX & new_mod);

/*------------------------------------------------------------*/
/* applies the above with all entries of new_mod              */
/*------------------------------------------------------------*/
void multimod(Vec<zz_pEX_augmented> & reduced, const zz_pEX_augmented & F, const zz_pX_CRT & new_mod);

/*------------------------------------------------------------*/
/* inverse of the previous one                                */
/* the moduli are the T1's of the entries of rems             */
/*------------------------------------------------------------*/
void combine(zz_pEX_augmented & result, const Vec<zz_pEX_augmented> & rems);

/*------------------------------------------------------------*/
/* splits F into polys where T2 is monic modulo T1            */
/*------------------------------------------------------------*/
void make_monic_and_split(Vec<zz_pEX_augmented> & monic_F, const zz_pEX_augmented & F);

/*------------------------------------------------------------*/
/* computes the squarefree part of f                          */
/* splits if needed, to ensure that the result is monic       */
/*------------------------------------------------------------*/
void squarefree(Vec<zz_pEX_augmented> & sqfree, const zz_pEX_augmented & f);

/*------------------------------------------------------------*/
/* computes the squarefree part of all entries of f           */
/* splits if needed, to ensure that the result is monic       */
/*------------------------------------------------------------*/
void squarefree(Vec<zz_pEX_augmented> & sqfree, const Vec<zz_pEX_augmented> & f);

/*------------------------------------------------------------*/
/* computes a decomposition of the zeros of f, g              */
/* assumes their resultant is non-zero                        */
/*------------------------------------------------------------*/
void solve(Vec<zz_pEX_augmented> & sols, const zz_pXY & f, const zz_pXY & g);

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* similar to zz_pEX_augmented, but with two polynomials in Y */
/* minimal functionalities are provided                       */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class pair_zz_pEX_augmented{
 public:
  zz_pX T1X; // keep a copy; cannot access it from T1 easily without switching moduli?
  zz_pEContext T1;
  zz_pEX T20, T21;

  pair_zz_pEX_augmented(){}

  pair_zz_pEX_augmented(const zz_pEX & t2_0, const zz_pEX & t2_1){
    T1X = zz_pE::modulus();
    T1 = zz_pEContext(T1X);
    T20 = t2_0;
    T21 = t2_1;
  }

  pair_zz_pEX_augmented(const zz_pX & t1, const zz_pXY & t2_0, const zz_pXY & t2_1){
    T1X = t1;
    T1 = zz_pEContext(t1);
    zz_pEPush push(T1);
    conv(T20, t2_0);
    conv(T21, t2_1);
  }

  pair_zz_pEX_augmented(const zz_pX & t1, const Vec<zz_pX> & t2_0, const Vec<zz_pX> & t2_1){
    T1X = t1;
    T1 = zz_pEContext(t1);
    zz_pEPush push(T1);
    Vec<zz_pE> t2E0, t2E1;

    conv(t2E0, t2_0);
    conv(T20, t2E0);
    conv(t2E1, t2_1);
    conv(T21, t2E1);
  }

  ~pair_zz_pEX_augmented(){}
};

/*------------------------------------------------------------*/
/* outputs as a sequence of M<y,x>                            */
/*------------------------------------------------------------*/
void magma_output(const pair_zz_pEX_augmented & F);

/*------------------------------------------------------------*/
/* assigns to "name" as a sequence of M<y,x>                  */
/*------------------------------------------------------------*/
void magma_assign(const pair_zz_pEX_augmented & F, const string & name);

/*------------------------------------------------------------*/
/* random, with deg(T1,x)=dX, deg(T2,y)=dY0, deg(T21,y)=dY1   */
/*------------------------------------------------------------*/
void random_pair_monic_zz_pEX_augmented(pair_zz_pEX_augmented & res, long dX, long dY0, long dY1);

/*------------------------------------------------------------*/
/* replaces T1 by new_mod and reduces T20 and T21 by it       */
/*------------------------------------------------------------*/
void mod(pair_zz_pEX_augmented & reduced, const pair_zz_pEX_augmented & orig, const zz_pX & new_mod);

/*------------------------------------------------------------*/
/* applies the above with all entries of new_mod              */
/*------------------------------------------------------------*/
void multimod(Vec<pair_zz_pEX_augmented> & reduced, const pair_zz_pEX_augmented & orig, const zz_pX_CRT & new_mod);

/*------------------------------------------------------------*/
/* inverse of the previous one                                */
/*------------------------------------------------------------*/
void combine(pair_zz_pEX_augmented & result, const Vec<pair_zz_pEX_augmented> & rems);

/*------------------------------------------------------------*/
/* splits F into polys where T2_idx is monic modulo T1        */
/*------------------------------------------------------------*/
void make_monic_and_split(Vec<pair_zz_pEX_augmented> & monic_F, const pair_zz_pEX_augmented & F, long idx);

/*------------------------------------------------------------*/
/* computes a regular gcd, splitting if needed                */
/*------------------------------------------------------------*/
void regular_gcd(Vec<zz_pEX_augmented> & h, const pair_zz_pEX_augmented & a_b);

#endif
