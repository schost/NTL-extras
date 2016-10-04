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
  Vec<BivariateModularComp> vec_M;
  cauchy_like_geometric_special CL,invA; //the inverse mod p
  Vec<ZZ_p> e, f; // for precondition the Cauchy Matrix
  ZZ c, d; // constant for the preconditioners of M
  Vec<ZZ_pX_Multipoint_FFT> vec_X_int, vec_Y_int; // to store already calculated powers
  Vec<ZZ> p_powers; // already calculated values of p^2^n
  BivariateModularComp *M; // the current M
  ZZ_pX_Multipoint_FFT *X_int, *Y_int; // the current preconditioner
  ZZX f_full_prec; // the polynomial in full precision
  Vec<long> type; // the type of the approximant
  long rank; // the rank of M
  long w; // w mod p
  long order; // order of w
  long p, prec, sizeX, sizeY;

  /**** helpers *****************************************/ 
  // switches the field to mod p^2^n
  void switch_context(long n);
  
  // calculates the order of w
  long find_order(zz_p w);

  // checks if every entry can be reconstructed in the current field
  bool can_reconstruct(const Vec<ZZ_p> &v, long n);

  // computes (D_e X_int) M (D_f Y_int)^t * b
  Vec<ZZ> mulA_right (Vec<ZZ_p> b);

  // computes x = Diagonal_d * a where Diagonal_d is the diagonal matrix constructed from d
  void mul_diagonal_right (Vec<ZZ_p> &x,const Vec<ZZ_p> &d, const Vec<ZZ_p> &a);
  
  // returns each block of v (based on the type) flipped
  Vec<ZZ_p> flip_on_type (const Vec<ZZ_p> &v);
  
  Vec<ZZ> find_original_sol(const Vec<ZZ_p> &b);
  
  // solves for Mx = b mod p^2^n
  void DAC (Vec<ZZ_p> &x, const Vec<ZZ_p> &b, long n);
  
  void reconstruct(Vec<Vec<ZZ>> &sol, const Vec<ZZ_p> &v, long n);
  
  Vec<Vec<ZZ>> flip_on_type(const Vec<Vec<ZZ>> &v);
  
  void dostuff();

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
