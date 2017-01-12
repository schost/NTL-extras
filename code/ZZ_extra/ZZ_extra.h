#ifndef __ZZ_EXTRA__H
#define __ZZ_EXTRA__H

#include <NTL/ZZ.h>
#include <gmp.h>

NTL_CLIENT

/*------------------------------------------------------------*/
/* ALLOC, SIZE, DATA = basic info on an _ntl_gbigint          */
/*------------------------------------------------------------*/
inline long ALLOC(_ntl_gbigint p)  { 
  return (((long *) p)[0]); 
}

inline long SIZE(_ntl_gbigint p) { 
  return (((long *) p)[1]); 
}

inline mp_limb_t * DATA(_ntl_gbigint p) { 
  return ((mp_limb_t *) (((long *) (p)) + 2)); 
}

inline long alloc(const ZZ& a)  {
  if (a == 0)
    return 0;
  return ALLOC(a.rep.rep);
}

inline long size(const ZZ& a) {
  if (a == 0)
    return 0;
  return SIZE(a.rep.rep);
}

inline mp_limb_t * data(const ZZ& a) {
  if (a == 0)
    return NULL;
  return DATA(a.rep.rep);
}

/*------------------------------------------------------------*/
/* multi-mod and CRT used for FFTs                            */
/*------------------------------------------------------------*/
void to_modular_rep(vec_long& x, const ZZ_p& a, const ZZ_pFFTInfoT *FFTInfo, ZZ_pTmpSpaceT *TmpSpace);
void from_modular_rep(ZZ_p& x, Vec<long>& avec, const ZZ_pFFTInfoT *FFTInfo, ZZ_pTmpSpaceT *TmpSpace);

#endif
