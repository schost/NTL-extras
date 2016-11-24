#ifndef __hermite_pade_h__
#define __hermite_pade_h__

#include <NTL/vector.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_pX.h>
#include <NTL/lzz_pX.h>
#include <NTL/tools.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include "mosaic_hankel.h"
#include "bivariatemodularcomp.h"
#include "ZZ_pX_CRT.h"
#include "cauchy_geometric_special.h"

// abstract base class for hermite pade classes
class hermite_pade{
	protected:
	NTL::Vec<long> type; // the type of the approximant
  long rank; // the rank of M
  long w; // w mod p
  long order; // order of w
  long p, prec, sizeX, sizeY;
  long added; // number of rows added
  long level;
  zz_pContext ctx; // zz_p 
  Vec<BivariateModularComp> vec_M;
  Vec<ZZ_pX_Multipoint_FFT> vec_X_int, vec_Y_int; // to store already calculated powers
  Vec<ZZ> p_powers; // already calculated values of p^2^n
  Vec<ZZ> e, f; // for precondition the Cauchy Matrix
  ZZ c, d; // constant for the preconditioners of M
  cauchy_like_geometric_special invA;
  
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
  
  // switches the field to mod p^2^n
  void switch_context(long n);
  
  virtual void set_up_bmc() = 0;

	void set_up_field(long fft_index);
	
	/*----------------------------------------------------------------*/
  /* multiplies b by the matrix CL =  D_e X_int M Y_int^t D_f       */
  /* (CL is cauchy-geometric-like)                                  */
  /* b need not have size CL.NumCols()                              */
  /*----------------------------------------------------------------*/
  Vec<ZZ_p> mulA_right (const Vec<ZZ_p> &b);
  
  /*----------------------------------------------------------------*/
  /* multiplies b by the matrix Y_int^t D_f                         */
  /*----------------------------------------------------------------*/
  Vec<ZZ_p> mul_Y_right(const Vec<ZZ_p> &b);
  
  /*----------------------------------------------------------------*/
  /* multiplies E_f * X_int * b                                     */
  /*----------------------------------------------------------------*/
  Vec<ZZ_p> mul_X_right(Vec<ZZ_p> b);

  /*----------------------------------------------------------------*/
  /* if Mx = b mod p^(2^{n-1}), updates x so that Mx = b mod p^(2^n)*/
  /*----------------------------------------------------------------*/
  void update_solution(Vec<ZZ>& x, const Vec<ZZ_p> &b, long n);
  
  /*----------------------------------------------------------------*/
  /* solves for Mx = b mod p^(2^n)                                  */
  /*----------------------------------------------------------------*/
  void DAC (Vec<ZZ> &x, const Vec<ZZ> &b, long n);
  
  /*----------------------------------------------------------------*/
  /* checks if every entry can be reconstructed in the current field*/
  /*----------------------------------------------------------------*/
  bool can_reconstruct(const Vec<ZZ_p> &v, long n);
  
  void reconstruct(Vec<Vec<ZZ>> &sol, const Vec<ZZ_p> &v, long n);
  
  /*----------------------------------------------------------------*/
  /* multiplies M*b                                                 */
  /*----------------------------------------------------------------*/
  virtual Vec<ZZ_p> mul_M_right(const Vec<ZZ_p> &b);
  
  hermite_pade(long fft_index): level{0}{
  	set_up_field(fft_index);
  }

  public:
  void find_rand_sol (NTL::Vec<NTL::Vec<NTL::ZZ>> &sol);
  virtual ~hermite_pade();
  
};

#endif













