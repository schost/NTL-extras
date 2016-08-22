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
  virtual void evaluate(Vec<Vec<ZZ_p>>& val, const Vec<ZZ_pX>& f) const = 0;
  virtual void interpolate(ZZ_pX& f, const Vec<ZZ_p>& val) const = 0;
 
  /*------------------------------------------------------------*/
  /* number of points                                           */
  /*------------------------------------------------------------*/
  long length() const{
    return n;
  }

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

  ZZ_pX_Multipoint_FFT(long n){
    k = NextPowerOfTwo(n);
    if (n != (1L << k)){
      cerr << "Wrong length for ZZ_pX_Multipoint_FFT\n";
      exit(-1);
    }

    this->n = n;
    if (ZZ_pInfo->p_info == NULL){
      cerr << "Attempt to init a ZZ_pX_Multipoint_FFT without ZZ_p::FFTInit\n";
      exit(-1);
    }
  }

  /*------------------------------------------------------------*/
  /* basic operations                                           */
  /*------------------------------------------------------------*/
  void evaluate(Vec<ZZ_p>& val, const ZZ_pX& f) const;
  void evaluate(Vec<Vec<ZZ_p>>& val, const Vec<ZZ_pX>& f) const;
  void interpolate(ZZ_pX& f, const Vec<ZZ_p>& val) const;
  
  ~ZZ_pX_Multipoint_FFT(){};

 private:
  long k;
};

/*------------------------------------------------------------*/
/* multiplicative order of a                                  */
/* -1 if a = 0                                                */
/*------------------------------------------------------------*/
long order(const ZZ_p& a);

/*------------------------------------------------------------*/
/* finds an element of order at least ord                     */
/* assumes it exists, does not verify                         */
/*------------------------------------------------------------*/
void element_of_order(ZZ_p& a, long ord);



#endif
