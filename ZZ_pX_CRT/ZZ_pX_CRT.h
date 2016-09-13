#ifndef ZZ_PX_CRT__H
#define ZZ_PX_CRT__H

#include <NTL/ZZ_pX.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* classes for CRT and multipoint evaluation over ZZ_p        */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* multipoint evaluation over ZZ_p                            */
/* abstract class                                             */
/* TODO: thresholds                                           */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class ZZ_pX_Multipoint{
 public:

  /*------------------------------------------------------------*/
  /* basic operations                                           */
  /*------------------------------------------------------------*/
  virtual void evaluate(Vec<ZZ_p>& val, const ZZ_pX& f) const = 0;
  virtual void interpolate(ZZ_pX& f, const Vec<ZZ_p>& val) const = 0;
 
  /*------------------------------------------------------------*/
  /* number of points                                           */
  /*------------------------------------------------------------*/
  long length() const{
    return n;
  }

  /*------------------------------------------------------------*/
  /* points                                                     */
  /*------------------------------------------------------------*/

  long n;
};


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* multipoint evaluation over ZZ_p                            */
/* FFT points                                                 */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class ZZ_pX_Multipoint_FFT : public ZZ_pX_Multipoint{
 public:

  ZZ_pX_Multipoint_FFT(){}

  ZZ_pX_Multipoint_FFT(const ZZ_p & w, long n);

  /*------------------------------------------------------------*/
  /* basic operations                                           */
  /*------------------------------------------------------------*/
  void evaluate(Vec<ZZ_p>& val, const ZZ_pX& f) const;

  void interpolate(ZZ_pX& f, const Vec<ZZ_p>& val) const{
    LogicError("inverse TFT not implemented");
  }
  
  ~ZZ_pX_Multipoint_FFT(){};

 private:
  long k, max_n;
  Vec<Vec<ZZ_p>> powersW;
  Vec<long> rev;
};

/*------------------------------------------------------------*/
/* finds a root of unity of order s mod p                     */
/* assumes s is a power of 2                                  */
/*------------------------------------------------------------*/
long find_root_of_unity(long p, long s);

/*------------------------------------------------------------*/
/* suppose that omega is a root of unity of order s mod p     */  
/* lifts omega to a root of unity Omega mod p^k               */
/*------------------------------------------------------------*/
void lift_root_of_unity(ZZ& Omega, long omega, long s, long p, long k);

#endif
