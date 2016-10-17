#ifndef __hermite_pade_h__
#define __hermite_pade_h__

#include "bivariatemodularcomp.h"
#include <NTL/ZZ_pX.h>
#include <NTL/lzz_pX.h>
#include <NTL/tools.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include "mosaic_hankel.h"
#include "ZZ_pX_CRT.h"
#include "cauchy_geometric_special.h"

class hermite_pade{

  long level;
  zz_pContext ctx; // zz_p 
  Vec<BivariateModularComp> vec_M;
  Vec<ZZ_pX_Multipoint_FFT> vec_X_int, vec_Y_int; // to store already calculated powers
  Vec<ZZ> p_powers; // already calculated values of p^2^n
  cauchy_like_geometric_special invA; //the inverse mod p
  Vec<ZZ> e, f; // for precondition the Cauchy Matrix
  ZZ c, d; // constant for the preconditioners of M

  ZZX f_full_prec; // the polynomial in full precision
  Vec<long> type; // the type of the approximant
  long rank; // the rank of M
  long w; // w mod p
  long order; // order of w
  long p, prec, sizeX, sizeY;
  long added; // number of rows added
  Vec<long> diagonals1, diagonals2;

  /**** helpers *****************************************/ 
  // switches the field to mod p^2^n
  void switch_context(long n);
  

  /*----------------------------------------------------------------*/
  /* returns a block of bi-sub-diagonal hankel matrices             */
  /* size is determined by the type                                 */
  /*----------------------------------------------------------------*/
  Vec<hankel> create_random_hankel_on_type();
  
  /*----------------------------------------------------------------*/
  /* multiplies b by the matrix CL =  D_e X_int M Y_int^t D_f       */
  /* (CL is cauchy-geometric-like)                                  */
  /* b need not have size CL.NumCols()                              */
  /*----------------------------------------------------------------*/
  Vec<ZZ_p> mulA_right (const Vec<ZZ_p> &b);

  /*----------------------------------------------------------------*/
  /* applies a block reversal to v                                  */
  /* e.g., type = [1,2] v = [v1,v2,v3,v4,v5]                        */
  /* returns [v2,v1,v4,v3,v2] (blocks have length type[i]+1)        */                                     
  /*----------------------------------------------------------------*/
  Vec<ZZ_p> flip_on_type (const Vec<ZZ_p> &v);

  /*----------------------------------------------------------------*/
  /* applies a block reversal to v                                  */
  /* e.g., type = [1,2] v = [v1,v2,v3,v4,v5]                        */
  /* returns [v2,v1,v4,v3,v2] (blocks have length type[i]+1)        */                                     
  /*----------------------------------------------------------------*/
  Vec<Vec<ZZ>> flip_on_type(const Vec<Vec<ZZ>> &v);
  
  /*----------------------------------------------------------------*/
  /* multiplies b by the matrix Y_int^t D_f                         */
  /*----------------------------------------------------------------*/
  Vec<ZZ_p> find_original_sol(const Vec<ZZ_p> &b);

  /*----------------------------------------------------------------*/
  /* if Mx = b mod p^(2^{n-1}), updates x so that Mx = b mod p^(2^n)*/
  /*----------------------------------------------------------------*/
  void update_solution(Vec<ZZ>& x, const Vec<ZZ_p> &b, long n);

  /*----------------------------------------------------------------*/
  /* computes Bv, where B is made of anti-diagonal matrices         */
  /*----------------------------------------------------------------*/
  Vec<ZZ_p> mul_bottom_mosaic_diagonal(const Vec<ZZ_p>& v);
  
  /*----------------------------------------------------------------*/
  /* solves for Mx = b mod p^(2^n)                                  */
  /*----------------------------------------------------------------*/
  void DAC (Vec<ZZ> &x, const Vec<ZZ> &b, long n);

  /*----------------------------------------------------------------*/
  /* checks if every entry can be reconstructed in the current field*/
  /*----------------------------------------------------------------*/
  bool can_reconstruct(const Vec<ZZ_p> &v, long n);
  
  void reconstruct(Vec<Vec<ZZ>> &sol, const Vec<ZZ_p> &v, long n);
 
  
  public:
  /** Ctor ***********************************************
  * f: the generating series                             *
  * type: the type of the approximant                    *
  *******************************************************/
  hermite_pade(const ZZX &f, const Vec<long>& type, long prec, long fft_init = 0);
  
  void find_rand_sol (Vec<Vec<ZZ>> &sol); 
 
  // computes a random solution to the system
  void get_rand_sol(Vec<ZZ> &);
};

#endif
