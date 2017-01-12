#ifndef ZZ_CRT__H
#define ZZ_CRT__H

#include <NTL/sp_arith.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include <gmp.h>

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* multiple remainder and CRT for integers                    */
/* (ripped off from NTL's code)                               */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* misc utility routines                                      */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* free _ntl_gbigints in a vector                             */
/*------------------------------------------------------------*/
void free_vec_gbigint(Vec<_ntl_gbigint>& v);

/*------------------------------------------------------------*/
/* set _ntl_gbigints to zero in a vector                      */
/*------------------------------------------------------------*/
void zero_vec_gbigint(Vec<_ntl_gbigint>& v);


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* classes for multi-mod                                      */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* base class                                                 */
/*------------------------------------------------------------*/
class ZZ_CRT_rem {
 public:
  long n;
  ZZ p;
  virtual void eval(Vec<long>& x_vec, const ZZ& a_ZZ) = 0;

 protected:

  Vec<long> moduli;
};

/*------------------------------------------------------------*/
/* ZZ_CRT_rem_tbl                                             */
/* DIRT: won't work if GMP has nails                          */
/* if defined(NTL_VIABLE_LL), defined(NTL_TBL_REM), (n <= 800)*/
/*------------------------------------------------------------*/
class ZZ_CRT_rem_tbl : public ZZ_CRT_rem {
 public:  
  ZZ_CRT_rem_tbl(){}
  ZZ_CRT_rem_tbl(const Vec<long>& moduli);
  void eval(Vec<long>& x_vec, const ZZ& a_ZZ);

 protected:
  Vec<mp_limb_t> inv_moduli;
  Mat<mp_limb_t> tbl;
};

/*------------------------------------------------------------*/
/* ZZ_CRT_rem_medium: class for large lengths                 */
/* (n >= 256), if rem_tbl not available                       */
/*------------------------------------------------------------*/
class ZZ_CRT_rem_fast : public ZZ_CRT_rem {
 public:
  ZZ_CRT_rem_fast(){}
  ZZ_CRT_rem_fast(const Vec<long>& moduli);
  void eval(Vec<long>& x_vec, const ZZ& a_ZZ);

  ~ZZ_CRT_rem_fast(){
    free_vec_gbigint(prod_vec);
    free_vec_gbigint(rem_vec);
  }
  
 protected:
  long levels;
  long vec_len;
  Vec<long> index_vec;
  Vec<_ntl_gbigint> prod_vec;
  Vec<_ntl_gbigint> rem_vec;
  long p_size;
};
/*------------------------------------------------------------*/
/* ZZ_CRT_rem_medium: class for medium lengths                */
/* (n >= 32 && n <= 256), if rem_tbl not available            */
/*------------------------------------------------------------*/
class ZZ_CRT_rem_medium : public ZZ_CRT_rem {
 public:
  ZZ_CRT_rem_medium(){}
  ZZ_CRT_rem_medium(const Vec<long>& moduli);
  void eval(Vec<long>& x_vec, const ZZ& a_ZZ);

  ~ZZ_CRT_rem_medium(){
    free_vec_gbigint(prod_vec);
    free_vec_gbigint(rem_vec);
  }

 protected:
  long levels;
  long vec_len;
  Vec<long> index_vec;
  Vec<long> len_vec, corr_vec;
  Vec<mulmod_precon_t> corraux_vec;
  Vec<mp_limb_t> inv_vec;
  Vec<_ntl_gbigint> prod_vec;
  Vec<_ntl_gbigint> rem_vec;
};

/*------------------------------------------------------------*/
/* ZZ_CRT_rem_medium: class for small lengths                 */
/* (n <= 1 && n < 32), if rem_tbl not available               */
/*------------------------------------------------------------*/
class ZZ_CRT_rem_basic : public ZZ_CRT_rem {
 public:
  ZZ_CRT_rem_basic(){}
  ZZ_CRT_rem_basic(const Vec<long>& moduli);
  void eval(Vec<long>& x_vec, const ZZ& a_ZZ);
};


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* classes for CRT                                            */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* base class                                                 */
/*------------------------------------------------------------*/
class ZZ_CRT_crt{
 public:
  virtual void crt(ZZ& a, const Vec<long>& b) = 0;
  
 protected:
  long n;
  ZZ p;
  Vec<long> moduli;
  Vec<double> moduli_recip;
  Vec<long> u;  // u[i] = (M/q[i])^{-1} mod q[i]
  Vec<mulmod_precon_t> uqinv;
  ZZ_ReduceStructAdapter reduce_struct;

  void make();
  void premultiply(Vec<long>& b, const Vec<long>& a);
  virtual void insert(long i, _ntl_gbigint m) = 0;
};

/*------------------------------------------------------------*/
/* ZZ_CRT_crt_tbl                                             */
/* if NTL_VIABLE_LL, NTL_CRT_ALTCODE, n <= 600                */
/* or NTL_VIABLE_LL, NTL_CRT_ALTCODE_SMALL, n <= 16           */
/*------------------------------------------------------------*/
class ZZ_CRT_crt_tbl : public ZZ_CRT_crt {
 public:
  ZZ_CRT_crt_tbl(){}
  ZZ_CRT_crt_tbl(const Vec<long>& moduli);
  void crt(ZZ& a, const Vec<long>& b);

 private:
  Vec<Vec<mp_limb_t>> v;
  long sz;
  void insert(long i, _ntl_gbigint m);
};

/*------------------------------------------------------------*/
/* ZZ_CRT_crt_fast                                            */
/* if n >= 600                                                */
/*------------------------------------------------------------*/
class ZZ_CRT_crt_fast : public ZZ_CRT_crt {
 public:
  ZZ_CRT_crt_fast(){}
  ZZ_CRT_crt_fast(const Vec<long>& moduli);
  void crt(ZZ& a, const Vec<long>& b);

  ~ZZ_CRT_crt_fast(){
    free_vec_gbigint(prod_vec);
    free_vec_gbigint(coeff_vec);
    free_vec_gbigint(rem_vec);
    free_vec_gbigint(temps);
  }    

 private:
   long levels;
   long vec_len;
   Vec<long> inv_vec;
   Vec<long> index_vec;
   Vec<_ntl_gbigint> prod_vec;
   Vec<_ntl_gbigint> coeff_vec;
   Vec<_ntl_gbigint> rem_vec;
   Vec<_ntl_gbigint> temps;
   Vec<long> val_vec;

   void insert(long i, _ntl_gbigint m);
};

/*------------------------------------------------------------*/
/* ZZ_CRT_crt_basic                                           */
/* all other cases                                            */
/*------------------------------------------------------------*/
class ZZ_CRT_crt_basic : public ZZ_CRT_crt {
 public:
  ZZ_CRT_crt_basic(){}
  ZZ_CRT_crt_basic(const Vec<long>& moduli);
  void crt(ZZ& a, const Vec<long>& b);

  ~ZZ_CRT_crt_basic(){
    free_vec_gbigint(v);
  }

 private:
   Vec<_ntl_gbigint> v;
   long sbuf;

   void insert(long i, _ntl_gbigint m);
};

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* a class that does everything                               */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

class ZZ_CRT {
 public:
  ZZ_CRT(const Vec<long>& moduli);

  void eval(Vec<long>& x, const ZZ& a){
    rem_struct->eval(x, a);
  }
  
  void crt(ZZ& a, const Vec<long>& x){
    crt_struct->crt(a, x);
  }

  // private:
  ZZ_CRT_crt* crt_struct;
  ZZ_CRT_rem* rem_struct;

  ZZ_CRT_rem_tbl rem_tbl;
  ZZ_CRT_rem_fast rem_fast;
  ZZ_CRT_rem_medium rem_medium;
  ZZ_CRT_rem_basic rem_basic;
  
#if ((defined(NTL_VIABLE_LL)) && (defined(NTL_CRT_ALTCODE))) ||  ((defined(NTL_VIABLE_LL)) && (defined(NTL_CRT_ALTCODE_SMALL)))
  ZZ_CRT_crt_tbl crt_tbl;
#endif
  ZZ_CRT_crt_fast crt_fast;
  ZZ_CRT_crt_basic crt_basic;

  /*------------------------------------------------------------*/
  /* helper routines to build attributes                        */
  /*------------------------------------------------------------*/
  void build_rem_struct(const Vec<long>& moduli);
  void build_crt_struct(const Vec<long>& moduli);
};


#endif
