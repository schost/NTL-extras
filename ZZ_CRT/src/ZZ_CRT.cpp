#include <NTL/sp_arith.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include <gmp.h>

#include "ZZ_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* multiple remainder and CRT for integers                    */
/* most of it is ripped off from NTL's code                   */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* some thresholds                                            */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
#define NTL_MAX_REM_TBL 800
#define NTL_MIN_REM_FAST 256
#define NTL_MIN_REM_MEDIUM 32

#define NTL_MAX_CRT_TBL 600

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* misc utility routines                                      */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* free _ntl_gbigints in a vector                             */
/*------------------------------------------------------------*/
void free_vec_gbigint(Vec<_ntl_gbigint>& v){
  for (long i = 0; i < v.length(); i++)
    _ntl_gfree(&v[i]);
}

/*------------------------------------------------------------*/
/* set _ntl_gbigints to zero in a vector                      */
/*------------------------------------------------------------*/
void zero_vec_gbigint(Vec<_ntl_gbigint>& v){
  for (long i = 0; i < v.length(); i++)
    v[i] = 0;
}

/*------------------------------------------------------------*/
/* ALLOC, SIZE, DATA = basic info on an _ntl_gbigint          */
/*------------------------------------------------------------*/
static
inline long& ALLOC(_ntl_gbigint p)  { 
  return (((long *) p)[0]); 
}

static
inline long& SIZE(_ntl_gbigint p) { 
  return (((long *) p)[1]); 
}

static
inline mp_limb_t * DATA(_ntl_gbigint p) { 
  return ((mp_limb_t *) (((long *) (p)) + 2)); 
}

/*------------------------------------------------------------*/
/* strips zeros, counts the size                              */
/*------------------------------------------------------------*/
static
inline void STRIP(long& sz, mp_limb_t *p) {
  long i;
  i = sz - 1;
  while (i >= 0 && p[i] == 0) i--;
  sz = i + 1;
}

/*------------------------------------------------------------*/
/* zero test                                                  */
/*------------------------------------------------------------*/
static
inline long ZEROP(_ntl_gbigint p) {
  return !p || !SIZE(p);
}

/*------------------------------------------------------------*/
/* decide whether new alloc is needed                         */
/*------------------------------------------------------------*/
static
inline long MustAlloc(_ntl_gbigint c, long len) { 
  return (!(c) || (ALLOC(c) >> 2) < (len)); 
}

/*------------------------------------------------------------*/
/* some kind of Newton iteration mod powers of 2?             */
/* DIRT: will not work with non-empty "nails"                 */
/*------------------------------------------------------------*/
mp_limb_t neg_inv_mod_limb(mp_limb_t m0) {
  mp_limb_t x; 
  long k;

  x = 1; 
  k = 1; 
  while (k < NTL_ZZ_NBITS) {
    x += x * (1 - x * m0);
    k <<= 1;
  }
  return - x;
}

/*------------------------------------------------------------*/
/* raises e to power p ?                                      */
/*------------------------------------------------------------*/
long SpecialPower(long e, long p) {
  long a;
  long x, y;

  a = (long) ((((mp_limb_t) 1) << (NTL_ZZ_NBITS-2)) % ((mp_limb_t) p));
  a = MulMod(a, 2, p);
  a = MulMod(a, 2, p);

  x = 1;
  y = a;
  while (e) {
    if (e & 1) x = MulMod(x, y, p);
    y = MulMod(y, y, p);
    e = e >> 1;
  }

  return x;
}

/*------------------------------------------------------------*/
/* small prime XGCD                                           */
/*------------------------------------------------------------*/
void sp_ext_eucl(long *dd, long *ss, long *tt, long a, long b){
   long  u, v, u0, v0, u1, v1, u2, v2, q, r;

   long aneg = 0, bneg = 0;

   if (a < 0) {
      if (a < -NTL_MAX_LONG) ResourceError("integer overflow");
      a = -a;
      aneg = 1;
   }

   if (b < 0) {
      if (b < -NTL_MAX_LONG) ResourceError("integer overflow");
      b = -b;
      bneg = 1;
   }

   u1=1; v1=0;
   u2=0; v2=1;
   u = a; v = b;

   while (v != 0) {
      q = u / v;
      r = u % v;
      u = v;
      v = r;
      u0 = u2;
      v0 = v2;
      u2 =  u1 - q*u2;
      v2 = v1- q*v2;
      u1 = u0;
      v1 = v0;
   }

   if (aneg)
      u1 = -u1;

   if (bneg)
      v1 = -v1;

   *dd = u;
   *ss = u1;
   *tt = v1;
}

/*------------------------------------------------------------*/
/* small prime modular inverse                                */
/*------------------------------------------------------------*/
long sp_inv_mod(long a, long n){
   long d, s, t;

   sp_ext_eucl(&d, &s, &t, a, n);
   if (d != 1) ArithmeticError("inverse undefined");
   if (s < 0)
      return s + n;
   else
      return s;
}

/*------------------------------------------------------------*/
/* Montgomery's redc                                          */
/*------------------------------------------------------------*/
void redc(_ntl_gbigint T, _ntl_gbigint N, long m, mp_limb_t inv, _ntl_gbigint res) {
  long n, sT, i;
  mp_limb_t *Ndata, *Tdata, *resdata, q, d, t, c;

  n = SIZE(N);
  Ndata = DATA(N);
  sT = SIZE(T);
  Tdata = DATA(T);
  resdata = DATA(res);

  for (i = sT; i < m+n; i++)
    Tdata[i] = 0;

  c = 0;
  for (i = 0; i < m; i++) {
    q = Tdata[i]*inv;
    d = mpn_addmul_1(Tdata+i, Ndata, n, q);

    t = Tdata[i+n] + d;
    Tdata[i+n] = t + c;
    if (t < d || (c == 1 && t + c  == 0)) 
      c = 1;
    else
      c = 0;
  }

  if (c) {
    mpn_sub_n(resdata, Tdata + m, Ndata, n);
  }
  else {
    for (i = 0; i < n; i++)
      resdata[i] = Tdata[m + i];
  }

  i = n;
  STRIP(i, resdata);

  SIZE(res) = i;
  SIZE(T) = 0;
}

/*------------------------------------------------------------*/
/* no idea what this does                                     */
/*------------------------------------------------------------*/
void gadd_mul_many(_ntl_gbigint *res, _ntl_gbigint *a, long *b, long n, long sz){
   mp_limb_t *xx, *yy; 
   long i, sx;
   long sy;
   mp_limb_t carry;

   sx = sz + 2;
   if (MustAlloc(*res, sx))
      _ntl_gsetlength(res, sx);

   xx = DATA(*res);

   for (i = 0; i < sx; i++)
      xx[i] = 0;

   for (i = 0; i < n; i++) {
      if (!a[i]) continue;

      yy = DATA(a[i]);
      sy = SIZE(a[i]); 

      if (!sy || !b[i]) continue;

      carry = mpn_addmul_1(xx, yy, sy, b[i]);
      yy = xx + sy;
      *yy += carry;

      if (*yy < carry) { /* unsigned comparison! */
         do {
            yy++;
            *yy += 1;
         } while (*yy == 0);
      }
   }

   while (sx > 0 && xx[sx-1] == 0) sx--;
   SIZE(*res) = sx;
}

/*------------------------------------------------------------*/
/* modular reduction 2 words -> 1 (probably)                  */
/*------------------------------------------------------------*/
static inline 
mp_limb_t tbl_red_21(mp_limb_t hi, mp_limb_t lo, long d, mp_limb_t dinv) {
  unsigned long H = (hi << (NTL_BITS_PER_LONG-NTL_SP_NBITS)) | (lo >> NTL_SP_NBITS);
  unsigned long Q = MulHiUL(H, dinv) + H;
  unsigned long rr = lo - Q*cast_unsigned(d); // rr in [0..4*d)
  long r = sp_CorrectExcess(rr, 2*d); // r in [0..2*d)
  r = sp_CorrectExcess(r, d);
  return r;
}

/*------------------------------------------------------------*/
/* modular reduction 3 words -> 1 (probably)                  */
/* NOTE: tbl_red_31 assumes x2 < d                            */
/*------------------------------------------------------------*/
static inline
mp_limb_t tbl_red_31(mp_limb_t x2, mp_limb_t x1, mp_limb_t x0, long d, mp_limb_t dinv) {
  mp_limb_t carry = tbl_red_21(x2, x1, d, dinv);
  return tbl_red_21(carry, x0, d, dinv);
}

/*------------------------------------------------------------*/
/* a simplified mod operation. Assumes                        */
/* a >= 0                                                     */
/* d > 0                                                      */
/* space for the result has already been allocated            */
/* inputs do not alias output                                 */
/*------------------------------------------------------------*/
void gmod_simple(_ntl_gbigint a, _ntl_gbigint d, _ntl_gbigint *rr) {
  _ntl_gbigint b = 0;

   long sa, sb, sd, sr;
   mp_limb_t *adata, *ddata, *bdata, *rdata;
   _ntl_gbigint r;

   if (ZEROP(a)) {
      _ntl_gzero(rr);
      return;
   }

   sa = SIZE(a);
   sd = SIZE(d);

   if (sa < sd) {
      _ntl_gcopy(a, rr);
      return;
   }

   sb = sa-sd+1;
   if (MustAlloc(b, sb))
      _ntl_gsetlength(&b, sb);

   sr = sd;
   r = *rr;

   adata = DATA(a);
   ddata = DATA(d);
   bdata = DATA(b);
   rdata = DATA(r);

   mpn_tdiv_qr(bdata, rdata, 0, adata, sa, ddata, sd);

   STRIP(sr, rdata);
   SIZE(r) = sr;

   _ntl_gfree(&b);
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* classes for multi-mod                                      */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* constructor of ZZ_CRT_rem_tbl                              */
/* does a lot of stuff I haven't read                         */
/*------------------------------------------------------------*/
ZZ_CRT_rem_tbl::ZZ_CRT_rem_tbl(const Vec<long>& q){

  if (q.length() > NTL_MAX_REM_TBL)
    return;
  
  long i, j;
  long qq, t, t1;

  moduli = q;
  n = moduli.length();

  p = to_ZZ(1);
  for (long i = 0; i < n; i++)
    p *= moduli[i];
  long sz = SIZE(p.rep);
   
  inv_moduli.SetLength(n);
  for (i = 0; i < n; i++) 
    inv_moduli[i] = (unsigned long) ( ((((NTL_ULL_TYPE) 1) << (NTL_SP_NBITS+NTL_BITS_PER_LONG))-1UL) / ((NTL_ULL_TYPE) moduli[i]) );
  
  tbl.SetDims(n, sz);
  
  for (i = 0; i < n; i++) {
    qq = moduli[i];
    t = 1;
    for (j = 0; j < NTL_ZZ_NBITS; j++) {
      t += t;
      if (t >= qq) t -= qq;
    }
    t1 = 1;
    tbl[i][0] = 1;
    for (j = 1; j < sz; j++) {
      t1 = MulMod(t1, t, qq);
      tbl[i][j] = t1;
    }
  }
}

/*------------------------------------------------------------*/
/* multiple remainder for ZZ_CRT_rem_tbl                      */
/* special case, some loop unrolling: slightly faster         */
/* DIRT: won't work if GMP has nails                          */
/*------------------------------------------------------------*/
#if (NTL_SP_NBITS == NTL_BITS_PER_LONG-2)
void ZZ_CRT_rem_tbl::eval(Vec<long>& x_vec, const ZZ& a_ZZ) {

  x_vec.SetLength(n);
  long * x = x_vec._vec__rep;
  ZZ a_rem = a_ZZ % p;
  _ntl_gbigint a = a_rem.rep;

  if (ZEROP(a)) {
    long i;
    for (i = 0; i < n; i++) x[i] = 0;
    return;
  }
  
  long sa = SIZE(a);
  mp_limb_t *adata = DATA(a);
  
  if (sa <= 4) {
    long i;
    for (i = 0; i < n; i++) {
      mp_limb_t *tp = tbl[i].rep; 
      NTL_ULL_TYPE acc = adata[0];
      long j;
      for (j = 1; j < sa; j++)
	acc += ((NTL_ULL_TYPE) adata[j]) * ((NTL_ULL_TYPE) tp[j]);
      
      mp_limb_t accvec[2];
      x[i] = tbl_red_31(0, acc >> NTL_BITS_PER_LONG, acc, moduli[i], inv_moduli[i]);
    }
  }
  else {
    long i;
    for (i = 0; i < n; i++) {
      mp_limb_t *ap = adata;
      mp_limb_t *tp = tbl[i].rep; 
      
      NTL_ULL_TYPE acc21;
      mp_limb_t acc0;
      
      {
	NTL_ULL_TYPE sum = ap[0];
	sum += ((NTL_ULL_TYPE) ap[1]) * ((NTL_ULL_TYPE) tp[1]);
	sum += ((NTL_ULL_TYPE) ap[2]) * ((NTL_ULL_TYPE) tp[2]);
	sum += ((NTL_ULL_TYPE) ap[3]) * ((NTL_ULL_TYPE) tp[3]);
	
	acc21 = sum >> NTL_BITS_PER_LONG;
	acc0 = sum;
      }
      
      long m=sa-4;
      ap += 4;
      tp += 4;
      
      for (; m >= 8; m -= 8, ap += 8, tp += 8) {
	{
	  NTL_ULL_TYPE sum = ((NTL_ULL_TYPE) ap[0]) * ((NTL_ULL_TYPE) tp[0]);
	  sum += ((NTL_ULL_TYPE) ap[1]) * ((NTL_ULL_TYPE) tp[1]);
	  sum += ((NTL_ULL_TYPE) ap[2]) * ((NTL_ULL_TYPE) tp[2]);
	  sum += ((NTL_ULL_TYPE) ap[3]) * ((NTL_ULL_TYPE) tp[3]);
	  
	  sum += acc0;
	  acc0 = sum;
	  acc21 += (mp_limb_t) (sum >> NTL_BITS_PER_LONG);
	}
	{
	  
	  NTL_ULL_TYPE sum = ((NTL_ULL_TYPE) ap[4+0]) * ((NTL_ULL_TYPE) tp[4+0]);
	  sum += ((NTL_ULL_TYPE) ap[4+1]) * ((NTL_ULL_TYPE) tp[4+1]);
	  sum += ((NTL_ULL_TYPE) ap[4+2]) * ((NTL_ULL_TYPE) tp[4+2]);
	  sum += ((NTL_ULL_TYPE) ap[4+3]) * ((NTL_ULL_TYPE) tp[4+3]);
	  
	  sum += acc0;
	  acc0 = sum;
	  acc21 += (mp_limb_t) (sum >> NTL_BITS_PER_LONG);
	}
      }
      
      for (; m >= 4; m -= 4, ap += 4, tp += 4) {
	NTL_ULL_TYPE sum = ((NTL_ULL_TYPE) ap[0]) * ((NTL_ULL_TYPE) tp[0]);
	sum += ((NTL_ULL_TYPE) ap[1]) * ((NTL_ULL_TYPE) tp[1]);
	sum += ((NTL_ULL_TYPE) ap[2]) * ((NTL_ULL_TYPE) tp[2]);
	sum += ((NTL_ULL_TYPE) ap[3]) * ((NTL_ULL_TYPE) tp[3]);
	
	sum += acc0;
	acc0 = sum;
	acc21 += (mp_limb_t) (sum >> NTL_BITS_PER_LONG);
      }
      
      if (m > 0) {
	NTL_ULL_TYPE sum = ((NTL_ULL_TYPE) ap[0]) * ((NTL_ULL_TYPE) tp[0]);
	for (m--, ap++, tp++; m > 0; m--, ap++, tp++)
	  sum += ((NTL_ULL_TYPE) ap[0]) * ((NTL_ULL_TYPE) tp[0]);
	
	sum += acc0;
	acc0 = sum;
	acc21 += (mp_limb_t) (sum >> NTL_BITS_PER_LONG);
      }
      
      x[i] = tbl_red_31(acc21 >> NTL_BITS_PER_LONG, acc21, acc0, moduli[i], inv_moduli[i]);
    }
  }
}

#else
/*----------------------------------------------------------------*/
/* General case: some loop unrolling (also using "Duff's Device") */
/* for the case where BPL-SPNBITS == 4: this is the common        */
/* case on 64-bit machines.  The loop unrolling and Duff seems    */
/* to shave off 5-10%                                             */
/* DIRT: won't work if GMP has nails                              */
/*----------------------------------------------------------------*/
#define TBL_UNROLL (1)
void ZZ_CRT_rem_tbl::eval(Vec<long>& x_vec, const ZZ& a_ZZ) {
  x_vec.SetLength(n);
  long * x = x_vec._vec__rep;
  ZZ a_rem = a_ZZ % p;
  _ntl_gbigint a = a_rem.rep;

  if (ZEROP(a)) {
    long i;
    for (i = 0; i < n; i++) x[i] = 0;
    return;
  }
  
  long sa = SIZE(a);
  mp_limb_t *adata = DATA(a);
  
  const long Bnd =  1L << (NTL_BITS_PER_LONG-NTL_SP_NBITS);
  
  if (sa <= Bnd) {
    long i;
    for (i = 0; i < n; i++) {
      mp_limb_t *tp = tbl[i]._vec__rep; 
      
      
      NTL_ULL_TYPE acc = adata[0];
      
#if (TBL_UNROLL && NTL_BITS_PER_LONG-NTL_SP_NBITS == 4)
      switch (sa) {
      case 16:  acc += ((NTL_ULL_TYPE) adata[16-1]) * ((NTL_ULL_TYPE) tp[16-1]);
      case 15:  acc += ((NTL_ULL_TYPE) adata[15-1]) * ((NTL_ULL_TYPE) tp[15-1]);
      case 14:  acc += ((NTL_ULL_TYPE) adata[14-1]) * ((NTL_ULL_TYPE) tp[14-1]);
      case 13:  acc += ((NTL_ULL_TYPE) adata[13-1]) * ((NTL_ULL_TYPE) tp[13-1]);
      case 12:  acc += ((NTL_ULL_TYPE) adata[12-1]) * ((NTL_ULL_TYPE) tp[12-1]);
      case 11:  acc += ((NTL_ULL_TYPE) adata[11-1]) * ((NTL_ULL_TYPE) tp[11-1]);
      case 10:  acc += ((NTL_ULL_TYPE) adata[10-1]) * ((NTL_ULL_TYPE) tp[10-1]);
      case 9:  acc += ((NTL_ULL_TYPE) adata[9-1]) * ((NTL_ULL_TYPE) tp[9-1]);
      case 8:  acc += ((NTL_ULL_TYPE) adata[8-1]) * ((NTL_ULL_TYPE) tp[8-1]);
      case 7:  acc += ((NTL_ULL_TYPE) adata[7-1]) * ((NTL_ULL_TYPE) tp[7-1]);
      case 6:  acc += ((NTL_ULL_TYPE) adata[6-1]) * ((NTL_ULL_TYPE) tp[6-1]);
      case 5:  acc += ((NTL_ULL_TYPE) adata[5-1]) * ((NTL_ULL_TYPE) tp[5-1]);
      case 4:  acc += ((NTL_ULL_TYPE) adata[4-1]) * ((NTL_ULL_TYPE) tp[4-1]);
      case 3:  acc += ((NTL_ULL_TYPE) adata[3-1]) * ((NTL_ULL_TYPE) tp[3-1]);
      case 2:  acc += ((NTL_ULL_TYPE) adata[2-1]) * ((NTL_ULL_TYPE) tp[2-1]);
      }
      
#else
      long j;
      for (j = 1; j < sa; j++)
	acc += ((NTL_ULL_TYPE) adata[j]) * ((NTL_ULL_TYPE) tp[j]);
#endif
      
      x[i] = tbl_red_31(0, acc >> NTL_ZZ_NBITS, acc, moduli[i], inv_moduli[i]);
    }
  }
  else {
    long i;
    for (i = 0; i < n; i++) {
      mp_limb_t *ap = adata;
      mp_limb_t *tp = tbl[i]._vec__rep; 
      
      NTL_ULL_TYPE acc21;
      mp_limb_t acc0;
      
      {
	NTL_ULL_TYPE sum = ap[0];
	
#if (TBL_UNROLL && NTL_BITS_PER_LONG-NTL_SP_NBITS == 4)
	sum += ((NTL_ULL_TYPE) ap[1]) * ((NTL_ULL_TYPE) tp[1]);
	sum += ((NTL_ULL_TYPE) ap[2]) * ((NTL_ULL_TYPE) tp[2]);
	sum += ((NTL_ULL_TYPE) ap[3]) * ((NTL_ULL_TYPE) tp[3]);
	sum += ((NTL_ULL_TYPE) ap[4]) * ((NTL_ULL_TYPE) tp[4]);
	sum += ((NTL_ULL_TYPE) ap[5]) * ((NTL_ULL_TYPE) tp[5]);
	sum += ((NTL_ULL_TYPE) ap[6]) * ((NTL_ULL_TYPE) tp[6]);
	sum += ((NTL_ULL_TYPE) ap[7]) * ((NTL_ULL_TYPE) tp[7]);
	sum += ((NTL_ULL_TYPE) ap[8]) * ((NTL_ULL_TYPE) tp[8]);
	sum += ((NTL_ULL_TYPE) ap[9]) * ((NTL_ULL_TYPE) tp[9]);
	sum += ((NTL_ULL_TYPE) ap[10]) * ((NTL_ULL_TYPE) tp[10]);
	sum += ((NTL_ULL_TYPE) ap[11]) * ((NTL_ULL_TYPE) tp[11]);
	sum += ((NTL_ULL_TYPE) ap[12]) * ((NTL_ULL_TYPE) tp[12]);
	sum += ((NTL_ULL_TYPE) ap[13]) * ((NTL_ULL_TYPE) tp[13]);
	sum += ((NTL_ULL_TYPE) ap[14]) * ((NTL_ULL_TYPE) tp[14]);
	sum += ((NTL_ULL_TYPE) ap[15]) * ((NTL_ULL_TYPE) tp[15]);
#else
	for (long j = 1; j < Bnd; j++)
	  sum += ((NTL_ULL_TYPE) ap[j]) * ((NTL_ULL_TYPE) tp[j]);
#endif
	
	acc21 = sum >> NTL_BITS_PER_LONG;
	acc0 = sum;
      }
      
      long m;
      for (m = sa-Bnd, ap += Bnd, tp += Bnd; m >= Bnd; m -= Bnd, ap += Bnd, tp += Bnd) {
	
	NTL_ULL_TYPE sum = ((NTL_ULL_TYPE) ap[0]) * ((NTL_ULL_TYPE) tp[0]);
	
#if (TBL_UNROLL && NTL_BITS_PER_LONG-NTL_SP_NBITS == 4)
	sum += ((NTL_ULL_TYPE) ap[1]) * ((NTL_ULL_TYPE) tp[1]);
	sum += ((NTL_ULL_TYPE) ap[2]) * ((NTL_ULL_TYPE) tp[2]);
	sum += ((NTL_ULL_TYPE) ap[3]) * ((NTL_ULL_TYPE) tp[3]);
	sum += ((NTL_ULL_TYPE) ap[4]) * ((NTL_ULL_TYPE) tp[4]);
	sum += ((NTL_ULL_TYPE) ap[5]) * ((NTL_ULL_TYPE) tp[5]);
	sum += ((NTL_ULL_TYPE) ap[6]) * ((NTL_ULL_TYPE) tp[6]);
	sum += ((NTL_ULL_TYPE) ap[7]) * ((NTL_ULL_TYPE) tp[7]);
	sum += ((NTL_ULL_TYPE) ap[8]) * ((NTL_ULL_TYPE) tp[8]);
	sum += ((NTL_ULL_TYPE) ap[9]) * ((NTL_ULL_TYPE) tp[9]);
	sum += ((NTL_ULL_TYPE) ap[10]) * ((NTL_ULL_TYPE) tp[10]);
	sum += ((NTL_ULL_TYPE) ap[11]) * ((NTL_ULL_TYPE) tp[11]);
	sum += ((NTL_ULL_TYPE) ap[12]) * ((NTL_ULL_TYPE) tp[12]);
	sum += ((NTL_ULL_TYPE) ap[13]) * ((NTL_ULL_TYPE) tp[13]);
	sum += ((NTL_ULL_TYPE) ap[14]) * ((NTL_ULL_TYPE) tp[14]);
	sum += ((NTL_ULL_TYPE) ap[15]) * ((NTL_ULL_TYPE) tp[15]);
#else
	for (long j = 1; j < Bnd; j++)
	  sum += ((NTL_ULL_TYPE) ap[j]) * ((NTL_ULL_TYPE) tp[j]);
#endif
	sum += acc0;
	acc0 = sum;
	acc21 += (mp_limb_t) (sum >> NTL_BITS_PER_LONG);
      }
      
      if (m > 0) {
	NTL_ULL_TYPE sum = ((NTL_ULL_TYPE) ap[0]) * ((NTL_ULL_TYPE) tp[0]);
	
#if (TBL_UNROLL && NTL_BITS_PER_LONG-NTL_SP_NBITS == 4)
	switch (m) {
	case 15:  sum += ((NTL_ULL_TYPE) ap[15-1]) * ((NTL_ULL_TYPE) tp[15-1]);
	case 14:  sum += ((NTL_ULL_TYPE) ap[14-1]) * ((NTL_ULL_TYPE) tp[14-1]);
	case 13:  sum += ((NTL_ULL_TYPE) ap[13-1]) * ((NTL_ULL_TYPE) tp[13-1]);
	case 12:  sum += ((NTL_ULL_TYPE) ap[12-1]) * ((NTL_ULL_TYPE) tp[12-1]);
	case 11:  sum += ((NTL_ULL_TYPE) ap[11-1]) * ((NTL_ULL_TYPE) tp[11-1]);
	case 10:  sum += ((NTL_ULL_TYPE) ap[10-1]) * ((NTL_ULL_TYPE) tp[10-1]);
	case 9:  sum += ((NTL_ULL_TYPE) ap[9-1]) * ((NTL_ULL_TYPE) tp[9-1]);
	case 8:  sum += ((NTL_ULL_TYPE) ap[8-1]) * ((NTL_ULL_TYPE) tp[8-1]);
	case 7:  sum += ((NTL_ULL_TYPE) ap[7-1]) * ((NTL_ULL_TYPE) tp[7-1]);
	case 6:  sum += ((NTL_ULL_TYPE) ap[6-1]) * ((NTL_ULL_TYPE) tp[6-1]);
	case 5:  sum += ((NTL_ULL_TYPE) ap[5-1]) * ((NTL_ULL_TYPE) tp[5-1]);
	case 4:  sum += ((NTL_ULL_TYPE) ap[4-1]) * ((NTL_ULL_TYPE) tp[4-1]);
	case 3:  sum += ((NTL_ULL_TYPE) ap[3-1]) * ((NTL_ULL_TYPE) tp[3-1]);
	case 2:  sum += ((NTL_ULL_TYPE) ap[2-1]) * ((NTL_ULL_TYPE) tp[2-1]);
	}
#else
	for (m--, ap++, tp++; m > 0; m--, ap++, tp++)
	  sum += ((NTL_ULL_TYPE) ap[0]) * ((NTL_ULL_TYPE) tp[0]);
#endif
	sum += acc0;
	acc0 = sum;
	acc21 += (mp_limb_t) (sum >> NTL_BITS_PER_LONG);
      }
      
      x[i] = tbl_red_31(acc21 >> NTL_BITS_PER_LONG, acc21, acc0, moduli[i], inv_moduli[i]);
    }
  }
}

#endif

/*------------------------------------------------------------*/
/* constructor of ZZ_CRT_rem_fast                             */
/* does a lot of stuff I haven't read                         */
/*------------------------------------------------------------*/
ZZ_CRT_rem_fast::ZZ_CRT_rem_fast(const Vec<long>& q){

  if (q.length() < NTL_MIN_REM_FAST)
    return;

  long i, j;
 
  moduli = q;
  n = moduli.length();

  p = to_ZZ(1);
  for (long i = 0; i < n; i++)
    p *= moduli[i];

  p_size = _ntl_gsize(p.rep);
  levels = 0;
  while ((n >> levels) >= 4) levels++;
  
  vec_len = (1L << levels) - 1;
  
  index_vec.SetLength(vec_len+1);
  prod_vec.SetLength(vec_len);
  for (long i = 0; i < vec_len; i++)
    prod_vec[i] = 0;
   
  index_vec[0] = 0;
  index_vec[1] = n;
  
  for (i = 0; i <= levels-2; i++) {
    long start = (1L << i) - 1;
    long finish = (1L << (i+1)) - 2;
    for (j = finish; j >= start; j--) {
      index_vec[2*j+2] = index_vec[j] + (index_vec[j+1] - index_vec[j])/2;
      index_vec[2*j+1] = index_vec[j];
    }
    index_vec[2*finish+3] = n;
  }

  for (i = (1L << (levels-1)) - 1; i < vec_len; i++) {
    /* multiply primes index_vec[i]..index_vec[i+1]-1 into 
     * prod_vec[i]
     */
    
    _ntl_gone(&prod_vec[i]);
    for (j = index_vec[i]; j < index_vec[i+1]; j++)
      _ntl_gsmul(prod_vec[i], moduli[j], &prod_vec[i]); 
  }
  
  for (i = (1L << (levels-1)) - 2; i >= 3; i--)
    _ntl_gmul(prod_vec[2*i+1], prod_vec[2*i+2], &prod_vec[i]);

  rem_vec.SetLength(vec_len);
  for (long i = 0; i < vec_len; i++)
    rem_vec[i] = 0;

  _ntl_gsetlength(&rem_vec[1], p_size);
  _ntl_gsetlength(&rem_vec[2], p_size);
  
  
  for (i = 1; i < (1L << (levels-1)) - 1; i++) {
    _ntl_gsetlength(&rem_vec[2*i+1], _ntl_gsize(prod_vec[2*i+1]));
    _ntl_gsetlength(&rem_vec[2*i+2], _ntl_gsize(prod_vec[2*i+2]));
   }
}

/*------------------------------------------------------------*/
/* multiple remainder for ZZ_CRT_rem_fast                     */
/*------------------------------------------------------------*/
void ZZ_CRT_rem_fast::eval(Vec<long>& x_vec, const ZZ& a_ZZ){

  x_vec.SetLength(n);
  long * x = x_vec._vec__rep;
  ZZ a_rem = a_ZZ % p;
  _ntl_gbigint a = a_rem.rep;

   long i, j;

   if (ZEROP(a)) {
      for (j = 0; j < n; j++)
         x[j] = 0;

      return;
   }

   _ntl_gcopy(a, &rem_vec[1]);
   _ntl_gcopy(a, &rem_vec[2]);

   for (i = 1; i < (1L << (levels-1)) - 1; i++) {
      gmod_simple(rem_vec[i], prod_vec[2*i+1], &rem_vec[2*i+1]);
      gmod_simple(rem_vec[i], prod_vec[2*i+2], &rem_vec[2*i+2]);
   }

   for (i = (1L << (levels-1)) - 1; i < vec_len; i++) {
      long lo = index_vec[i];
      long hi = index_vec[i+1];
      mp_limb_t *s1p = DATA(rem_vec[i]);

      long s1size = SIZE(rem_vec[i]);
      if (s1size == 0) {
         for (j = lo; j <hi; j++)
            x[j] = 0;
      }
      else {
         for (j = lo; j < hi; j++)
            x[j] = mpn_mod_1(s1p, s1size, moduli[j]);
      }
   }

}

/*------------------------------------------------------------*/
/* constructor of ZZ_CRT_rem_medium                           */
/* does a lot of stuff I haven't read                         */
/*------------------------------------------------------------*/
ZZ_CRT_rem_medium::ZZ_CRT_rem_medium(const Vec<long>& q){

  if (q.length() < NTL_MIN_REM_MEDIUM)
    return;
  if (q.length() >= NTL_MIN_REM_FAST)
    return;

  long i, j;
 
  moduli = q;
  n = moduli.length();

  p = to_ZZ(1);
  for (long i = 0; i < n; i++)
    p *= moduli[i];

  levels = 0;
  while ((n >> levels) >= 4) 
    levels++;
  vec_len = (1L << levels) - 1;

  index_vec.SetLength(vec_len+1);
  len_vec.SetLength(vec_len);
  inv_vec.SetLength(vec_len);
  corr_vec.SetLength(n);
  corraux_vec.SetLength(n);

  prod_vec.SetLength(vec_len);
  zero_vec_gbigint(prod_vec);
      
  index_vec[0] = 0;
  index_vec[1] = n;
  
  for (i = 0; i <= levels-2; i++) {
    long start = (1L << i) - 1;
    long finish = (1L << (i+1)) - 2;
    for (j = finish; j >= start; j--) {
      index_vec[2*j+2] = index_vec[j] + (index_vec[j+1] - index_vec[j])/2;
      index_vec[2*j+1] = index_vec[j];
    }
    index_vec[2*finish+3] = n;
  }

  for (i = (1L << (levels-1)) - 1; i < vec_len; i++) {
    /* multiply primes index_vec[i]..index_vec[i+1]-1 into 
     * prod_vec[i]
     */
    _ntl_gone(&prod_vec[i]);
    for (j = index_vec[i]; j < index_vec[i+1]; j++)
      _ntl_gsmul(prod_vec[i], moduli[j], &prod_vec[i]); 
  }

  for (i = (1L << (levels-1)) - 2; i >= 3; i--){
    _ntl_gmul(prod_vec[2*i+1], prod_vec[2*i+2], &prod_vec[i]);
  }

  for (i = 3; i < vec_len; i++)
    len_vec[i] = _ntl_gsize(prod_vec[i]);

  /* Set len_vec[1] = len_vec[2] = 
   *    max(_ntl_gsize(modulus), len_vec[3..6]).
   * This is a bit paranoid, but it makes the code
   * more robust. */
  
  j = _ntl_gsize(p.rep);
  for (i = 3; i <= 6; i++)
    if (len_vec[i] > j) j = len_vec[i];
  len_vec[1] = len_vec[2] = j;

  for (i = 3; i < vec_len; i++)
    inv_vec[i] = neg_inv_mod_limb(DATA(prod_vec[i])[0]);


  for (i = (1L << (levels-1)) - 1; i < vec_len; i++) {
    for (j = index_vec[i]; j < index_vec[i+1]; j++) {
      corr_vec[j] = SpecialPower(len_vec[1] - len_vec[i], moduli[j]);
      corraux_vec[j] = PrepMulModPrecon(corr_vec[j], moduli[j]);
    }
  }

  /* allocate length in advance to streamline eval code */
  rem_vec.SetLength(vec_len);
  for (long i = 0; i < vec_len; i++)
    rem_vec[i] = 0;

  _ntl_gsetlength(&rem_vec[0], len_vec[1]); /* a special temp */
  for (i = 1; i < vec_len; i++)
    _ntl_gsetlength(&rem_vec[i], len_vec[i]);
}

/*------------------------------------------------------------*/
/* multiple remainder for ZZ_CRT_rem_medium                   */
/*------------------------------------------------------------*/
void ZZ_CRT_rem_medium::eval(Vec<long>& x_vec, const ZZ& a_ZZ){

  x_vec.SetLength(n);
  long * x = x_vec._vec__rep;
  ZZ a_rem = a_ZZ % p;
  _ntl_gbigint a = a_rem.rep;

  long i, j;

  if (ZEROP(a)) {
    for (j = 0; j < n; j++)
      x[j] = 0;

    return;
  }

  _ntl_gcopy(a, &rem_vec[1]);
  _ntl_gcopy(a, &rem_vec[2]);

  for (i = 1; i < (1L << (levels-1)) - 1; i++) {
    _ntl_gcopy(rem_vec[i], &rem_vec[0]);
    redc(rem_vec[0], prod_vec[2*i+1], len_vec[i]-len_vec[2*i+1],
	 inv_vec[2*i+1], rem_vec[2*i+1]);
    redc(rem_vec[i], prod_vec[2*i+2], len_vec[i]-len_vec[2*i+2],
	 inv_vec[2*i+2], rem_vec[2*i+2]);
  }

  for (i = (1L << (levels-1)) - 1; i < vec_len; i++) {
    long lo = index_vec[i];
    long hi = index_vec[i+1];
    mp_limb_t *s1p = DATA(rem_vec[i]);
    long s1size = SIZE(rem_vec[i]);
    if (s1size == 0) {
      for (j = lo; j <hi; j++)
	x[j] = 0;
    }
    else {
      for (j = lo; j < hi; j++) {
	long t = mpn_mod_1(s1p, s1size, moduli[j]);
	x[j] = MulModPrecon(t, corr_vec[j], moduli[j], corraux_vec[j]);
      }
    }
  }
}

/*------------------------------------------------------------*/
/* constructor of ZZ_CRT_rem_basic                            */
/* does a lot of stuff I haven't read                         */
/*------------------------------------------------------------*/
ZZ_CRT_rem_basic::ZZ_CRT_rem_basic(const Vec<long>& q){
  
  if (q.length() >= NTL_MIN_REM_MEDIUM)
    return;

  moduli = q;
  n = moduli.length();
  
  p = to_ZZ(1);
  for (long i = 0; i < n; i++)
    p *= moduli[i];
}


/*------------------------------------------------------------*/
/* multiple remainder for ZZ_CRT_rem_basic                    */
/*------------------------------------------------------------*/
void ZZ_CRT_rem_basic::eval(Vec<long>& x_vec, const ZZ& a_ZZ){
  x_vec.SetLength(n);
  long * x = x_vec._vec__rep;
  ZZ a_rem = a_ZZ % p;
  _ntl_gbigint a = a_rem.rep;

 
  long j;
  mp_limb_t *adata;
  long sa;
  
  if (!a) 
    sa = 0;
  else
    sa = SIZE(a);
  
  if (sa == 0) {
    for (j = 0; j < n; j++)
      x[j] = 0;
    return;
  }
  
  adata = DATA(a);
  
  for (j = 0; j < n; j++)
    x[j] = mpn_mod_1(adata, sa, moduli[j]);
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* classes for CRT                                            */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* base class methods                                         */
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* does a few precomputations                                 */
/* useful only for basic and tbl                              */
/*------------------------------------------------------------*/
void ZZ_CRT_crt::make(){
  moduli_recip.SetLength(n);
  u.SetLength(n);
  uqinv.SetLength(n);
  // montgomery
  reduce_struct.init(p, ZZ(n) << NTL_SP_NBITS);
  
  for (long i = 0; i < n; i++) {
    long q = moduli[i];
    mulmod_t qinv = PrepMulMod(q);
    
    ZZ M2;
    div(M2, p, q);  // = (M/q) rem p
    long t = rem(M2, q);
    t = InvMod(t, q);
    
    // montgomery
    reduce_struct.adjust(M2);
    insert(i, M2.rep);
    
    moduli_recip[i] = 1/double(q);
    u[i] = t;
    uqinv[i] = PrepMulModPrecon(t, q, qinv);
  }
}

/*------------------------------------------------------------*/
/* a pre-multiplication step                                  */
/* useful only for basic and tbl                              */
/*------------------------------------------------------------*/
void ZZ_CRT_crt::premultiply(Vec<long>& b, const Vec<long>& a){
  b.SetLength(n);
  for (long i = 0; i < n; i++) {
    long r = MulModPrecon(a[i], u[i], moduli[i], uqinv[i]);
    b[i] = r;
  }
}

/*------------------------------------------------------------*/
/* constructor of ZZ_CRT_crt_tbl                              */
/* does a lot of stuff I haven't read                         */
/* call to make at the end (use Montgomery)                   */
/*------------------------------------------------------------*/
ZZ_CRT_crt_tbl::ZZ_CRT_crt_tbl(const Vec<long>& q){
  if (q.length() > NTL_MAX_CRT_TBL)
    return;

  moduli = q;
  n = moduli.length();

  p = to_ZZ(1);
  for (long i = 0; i < n; i++)
    p *= moduli[i];
  
  sz = SIZE(p.rep);
  v.SetLength(sz);
  for (long i = 0; i < sz; i++)
    v[i].SetLength(n);

  make(); 
}

/*------------------------------------------------------------*/
/* inserts m into v[i]                                        */
/*------------------------------------------------------------*/
void ZZ_CRT_crt_tbl::insert(long i, _ntl_gbigint m){
  if (!m) 
    for (long j = 0; j < sz; j++) v[j][i] = 0;
  else {
    long sm = SIZE(m);
    if (sm < 0 || sm > sz) 
      LogicError("insert: bad args");
    const mp_limb_t *mdata = DATA(m);
    for (long j = 0; j < sm; j++) 
      v[j][i] = mdata[j];
    for (long j = sm; j < sz; j++)
      v[j][i] = 0;
  }
}

/*------------------------------------------------------------*/
/* CRT for ZZ_CRT_crt_tbl                                     */
/* premultiply at the beginning, reduce.eval at the end       */
/*------------------------------------------------------------*/
#define CRT_ALTCODE_UNROLL (1)
void ZZ_CRT_crt_tbl::crt(ZZ& x_ZZ, const Vec<long>& bb_in){

  Vec<long> b_in;
  premultiply(b_in, bb_in);

  ZZ xx1;
  _ntl_gbigint *x = &xx1.rep;
  long sx;
  _ntl_gbigint x1;
  long i;
  
  long* b = b_in._vec__rep;
  
  // quick test for zero vector
  // most likely, they are either all zero (if we are working 
  // with some sparse polynomials) or none of them are zero,
  // so in the general case, this should go fast
  if (!b_in[0]) {
    i = 1;
    while (i < n && !b_in[i]) i++;
    if (i >= n) {
      x_ZZ = to_ZZ(0);
      return;
    }
  }

  sx = sz + 2;
  _ntl_gsetlength(x, sx);
  x1 = *x;
  mp_limb_t * NTL_RESTRICT xx = DATA(x1);
  
  
  const long Bnd = 1L << (NTL_BITS_PER_LONG-NTL_SP_NBITS);
  
  if (n <= Bnd) {
    mp_limb_t carry=0;
    
    for (i = 0; i < sz; i++) {
      const mp_limb_t *row = v[i]._vec__rep;
      
      NTL_ULL_TYPE acc = ((NTL_ULL_TYPE) row[0]) * ((NTL_ULL_TYPE) (mp_limb_t) b[0]);
      
#if (CRT_ALTCODE_UNROLL && NTL_BITS_PER_LONG-NTL_SP_NBITS == 4)
      switch (n) {
      case 16: acc += ((NTL_ULL_TYPE) row[16-1]) * ((NTL_ULL_TYPE) (mp_limb_t) b[16-1]);
      case 15: acc += ((NTL_ULL_TYPE) row[15-1]) * ((NTL_ULL_TYPE) (mp_limb_t) b[15-1]);
      case 14: acc += ((NTL_ULL_TYPE) row[14-1]) * ((NTL_ULL_TYPE) (mp_limb_t) b[14-1]);
      case 13: acc += ((NTL_ULL_TYPE) row[13-1]) * ((NTL_ULL_TYPE) (mp_limb_t) b[13-1]);
      case 12: acc += ((NTL_ULL_TYPE) row[12-1]) * ((NTL_ULL_TYPE) (mp_limb_t) b[12-1]);
      case 11: acc += ((NTL_ULL_TYPE) row[11-1]) * ((NTL_ULL_TYPE) (mp_limb_t) b[11-1]);
      case 10: acc += ((NTL_ULL_TYPE) row[10-1]) * ((NTL_ULL_TYPE) (mp_limb_t) b[10-1]);
      case 9: acc += ((NTL_ULL_TYPE) row[9-1]) * ((NTL_ULL_TYPE) (mp_limb_t) b[9-1]);
      case 8: acc += ((NTL_ULL_TYPE) row[8-1]) * ((NTL_ULL_TYPE) (mp_limb_t) b[8-1]);
      case 7: acc += ((NTL_ULL_TYPE) row[7-1]) * ((NTL_ULL_TYPE) (mp_limb_t) b[7-1]);
      case 6: acc += ((NTL_ULL_TYPE) row[6-1]) * ((NTL_ULL_TYPE) (mp_limb_t) b[6-1]);
      case 5: acc += ((NTL_ULL_TYPE) row[5-1]) * ((NTL_ULL_TYPE) (mp_limb_t) b[5-1]);
      case 4: acc += ((NTL_ULL_TYPE) row[4-1]) * ((NTL_ULL_TYPE) (mp_limb_t) b[4-1]);
      case 3: acc += ((NTL_ULL_TYPE) row[3-1]) * ((NTL_ULL_TYPE) (mp_limb_t) b[3-1]);
      case 2: acc += ((NTL_ULL_TYPE) row[2-1]) * ((NTL_ULL_TYPE) (mp_limb_t) b[2-1]);
      }
#else
      for (j = 1; j < n; j++) 
	acc += ((NTL_ULL_TYPE) row[j]) * ((NTL_ULL_TYPE) (mp_limb_t) b[j]);
#endif
      
      acc += carry;
      xx[i] = acc;
      carry = acc >> NTL_BITS_PER_LONG;
    }
    
    xx[sz] = carry;
    xx[sz+1] = 0;
  }
  else {
    NTL_ULL_TYPE carry=0;
    
    for (i = 0; i < sz; i++) {
      const mp_limb_t *row = v[i]._vec__rep;
      
      NTL_ULL_TYPE acc21;
      mp_limb_t acc0;
      
      {
	NTL_ULL_TYPE sum = ((NTL_ULL_TYPE) row[0]) * ((NTL_ULL_TYPE) (mp_limb_t) b[0]);
	
#if (CRT_ALTCODE_UNROLL && NTL_BITS_PER_LONG-NTL_SP_NBITS == 4)
	sum += ((NTL_ULL_TYPE) row[1]) * ((NTL_ULL_TYPE) (mp_limb_t) b[1]);
	sum += ((NTL_ULL_TYPE) row[2]) * ((NTL_ULL_TYPE) (mp_limb_t) b[2]);
	sum += ((NTL_ULL_TYPE) row[3]) * ((NTL_ULL_TYPE) (mp_limb_t) b[3]);
	sum += ((NTL_ULL_TYPE) row[4]) * ((NTL_ULL_TYPE) (mp_limb_t) b[4]);
	sum += ((NTL_ULL_TYPE) row[5]) * ((NTL_ULL_TYPE) (mp_limb_t) b[5]);
	sum += ((NTL_ULL_TYPE) row[6]) * ((NTL_ULL_TYPE) (mp_limb_t) b[6]);
	sum += ((NTL_ULL_TYPE) row[7]) * ((NTL_ULL_TYPE) (mp_limb_t) b[7]);
	sum += ((NTL_ULL_TYPE) row[8]) * ((NTL_ULL_TYPE) (mp_limb_t) b[8]);
	sum += ((NTL_ULL_TYPE) row[9]) * ((NTL_ULL_TYPE) (mp_limb_t) b[9]);
	sum += ((NTL_ULL_TYPE) row[10]) * ((NTL_ULL_TYPE) (mp_limb_t) b[10]);
	sum += ((NTL_ULL_TYPE) row[11]) * ((NTL_ULL_TYPE) (mp_limb_t) b[11]);
	sum += ((NTL_ULL_TYPE) row[12]) * ((NTL_ULL_TYPE) (mp_limb_t) b[12]);
	sum += ((NTL_ULL_TYPE) row[13]) * ((NTL_ULL_TYPE) (mp_limb_t) b[13]);
	sum += ((NTL_ULL_TYPE) row[14]) * ((NTL_ULL_TYPE) (mp_limb_t) b[14]);
	sum += ((NTL_ULL_TYPE) row[15]) * ((NTL_ULL_TYPE) (mp_limb_t) b[15]);
#elif (CRT_ALTCODE_UNROLL && NTL_BITS_PER_LONG-NTL_SP_NBITS == 2)
	sum += ((NTL_ULL_TYPE) row[1]) * ((NTL_ULL_TYPE) (mp_limb_t) b[1]);
	sum += ((NTL_ULL_TYPE) row[2]) * ((NTL_ULL_TYPE) (mp_limb_t) b[2]);
	sum += ((NTL_ULL_TYPE) row[3]) * ((NTL_ULL_TYPE) (mp_limb_t) b[3]);
#else
	for (j = 1; j < Bnd; j++)
	  sum += ((NTL_ULL_TYPE) row[j]) * ((NTL_ULL_TYPE) (mp_limb_t) b[j]);
#endif
	
	acc21 = sum >> NTL_BITS_PER_LONG;
	acc0 = sum;
      }
      
      const mp_limb_t *ap = row;
      const long *tp = b;
      
#if (CRT_ALTCODE_UNROLL && NTL_BITS_PER_LONG-NTL_SP_NBITS == 2)
      long m = n - 4;
      ap += 4;
      tp += 4;
      
      for (; m >= 8; m -= 8, ap += 8, tp += 8) {
	{
	  NTL_ULL_TYPE sum = ((NTL_ULL_TYPE) ap[0]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[0]);
	  sum += ((NTL_ULL_TYPE) ap[1]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[1]);
	  sum += ((NTL_ULL_TYPE) ap[2]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[2]);
	  sum += ((NTL_ULL_TYPE) ap[3]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[3]);
	  sum += acc0;
	  acc0 = sum;
	  acc21 += (mp_limb_t) (sum >> NTL_BITS_PER_LONG);
	}
	{
	  NTL_ULL_TYPE sum = ((NTL_ULL_TYPE) ap[4+0]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[4+0]);
	  sum += ((NTL_ULL_TYPE) ap[4+1]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[4+1]);
	  sum += ((NTL_ULL_TYPE) ap[4+2]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[4+2]);
	  sum += ((NTL_ULL_TYPE) ap[4+3]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[4+3]);
	  sum += acc0;
	  acc0 = sum;
	  acc21 += (mp_limb_t) (sum >> NTL_BITS_PER_LONG);
	}
      }
      
      for (; m >= 4; m -= 4, ap += 4, tp += 4) {
	NTL_ULL_TYPE sum = ((NTL_ULL_TYPE) ap[0]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[0]);
	sum += ((NTL_ULL_TYPE) ap[1]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[1]);
	sum += ((NTL_ULL_TYPE) ap[2]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[2]);
	sum += ((NTL_ULL_TYPE) ap[3]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[3]);
	sum += acc0;
	acc0 = sum;
	acc21 += (mp_limb_t) (sum >> NTL_BITS_PER_LONG);
      }
      
      
#else
      long m;
      for (m = n-Bnd, ap += Bnd, tp += Bnd; m >= Bnd; m -= Bnd, ap += Bnd, tp += Bnd) {
	
	NTL_ULL_TYPE sum = ((NTL_ULL_TYPE) ap[0]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[0]);
	
#if (CRT_ALTCODE_UNROLL && NTL_BITS_PER_LONG-NTL_SP_NBITS == 4)
	sum += ((NTL_ULL_TYPE) ap[1]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[1]);
	sum += ((NTL_ULL_TYPE) ap[2]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[2]);
	sum += ((NTL_ULL_TYPE) ap[3]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[3]);
	sum += ((NTL_ULL_TYPE) ap[4]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[4]);
	sum += ((NTL_ULL_TYPE) ap[5]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[5]);
	sum += ((NTL_ULL_TYPE) ap[6]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[6]);
	sum += ((NTL_ULL_TYPE) ap[7]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[7]);
	sum += ((NTL_ULL_TYPE) ap[8]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[8]);
	sum += ((NTL_ULL_TYPE) ap[9]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[9]);
	sum += ((NTL_ULL_TYPE) ap[10]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[10]);
	sum += ((NTL_ULL_TYPE) ap[11]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[11]);
	sum += ((NTL_ULL_TYPE) ap[12]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[12]);
	sum += ((NTL_ULL_TYPE) ap[13]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[13]);
	sum += ((NTL_ULL_TYPE) ap[14]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[14]);
	sum += ((NTL_ULL_TYPE) ap[15]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[15]);
#else
	for (long j = 1; j < Bnd; j++)
	  sum += ((NTL_ULL_TYPE) ap[j]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[j]);
#endif
	
	sum += acc0;
	acc0 = sum;
	acc21 += (mp_limb_t) (sum >> NTL_BITS_PER_LONG);
      }
#endif
      
      if (m > 0) {
	NTL_ULL_TYPE sum = ((NTL_ULL_TYPE) ap[0]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[0]);
	
#if (CRT_ALTCODE_UNROLL && NTL_BITS_PER_LONG-NTL_SP_NBITS == 4)
	switch (m) {
	case 15:  sum += ((NTL_ULL_TYPE) ap[15-1]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[15-1]);
	case 14:  sum += ((NTL_ULL_TYPE) ap[14-1]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[14-1]);
	case 13:  sum += ((NTL_ULL_TYPE) ap[13-1]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[13-1]);
	case 12:  sum += ((NTL_ULL_TYPE) ap[12-1]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[12-1]);
	case 11:  sum += ((NTL_ULL_TYPE) ap[11-1]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[11-1]);
	case 10:  sum += ((NTL_ULL_TYPE) ap[10-1]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[10-1]);
	case 9:  sum += ((NTL_ULL_TYPE) ap[9-1]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[9-1]);
	case 8:  sum += ((NTL_ULL_TYPE) ap[8-1]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[8-1]);
	case 7:  sum += ((NTL_ULL_TYPE) ap[7-1]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[7-1]);
	case 6:  sum += ((NTL_ULL_TYPE) ap[6-1]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[6-1]);
	case 5:  sum += ((NTL_ULL_TYPE) ap[5-1]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[5-1]);
	case 4:  sum += ((NTL_ULL_TYPE) ap[4-1]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[4-1]);
	case 3:  sum += ((NTL_ULL_TYPE) ap[3-1]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[3-1]);
	case 2:  sum += ((NTL_ULL_TYPE) ap[2-1]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[2-1]);
	}
#else
	for (m--, ap++, tp++; m > 0; m--, ap++, tp++)
	  sum += ((NTL_ULL_TYPE) ap[0]) * ((NTL_ULL_TYPE) (mp_limb_t) tp[0]);
#endif
	sum += acc0;
	acc0 = sum;
	acc21 += (mp_limb_t) (sum >> NTL_BITS_PER_LONG);
      }
      
      carry += acc0;
      xx[i] = carry;
      acc21 += ((mp_limb_t) (carry >> NTL_BITS_PER_LONG));
      carry = acc21;
    }
    
    xx[sz] = carry;
    xx[sz+1] = carry >> NTL_BITS_PER_LONG;
  }
    
  while (sx > 0 && xx[sx-1] == 0) 
    sx--;
  SIZE(x1) = sx;

  reduce_struct.eval(x_ZZ, xx1);
}


/*------------------------------------------------------------*/
/* constructor of ZZ_CRT_crt_fast                             */
/* does a lot of stuff I haven't read                         */
/*------------------------------------------------------------*/
#define GCRT_TMPS (2)
ZZ_CRT_crt_fast::ZZ_CRT_crt_fast(const Vec<long>& q){
  if (q.length() < NTL_MAX_CRT_TBL)
    return;
  
  moduli = q;
  n = moduli.length();

  p = to_ZZ(1);
  for (long i = 0; i < n; i++)
    p *= moduli[i];
  
  levels = 0;
  while ((n >> levels) >= 16) 
    levels++;
  vec_len = (1L << levels) - 1;

  val_vec.SetLength(n);
  temps.SetLength(GCRT_TMPS);
  zero_vec_gbigint(temps);
  rem_vec.SetLength(vec_len);
  zero_vec_gbigint(rem_vec);
  inv_vec.SetLength(n);
  index_vec.SetLength(vec_len+1);
  prod_vec.SetLength(vec_len);
  zero_vec_gbigint(prod_vec);
  coeff_vec.SetLength(n);
  zero_vec_gbigint(coeff_vec);
  
  index_vec[0] = 0;
  index_vec[1] = n;

  for (long i = 0; i <= levels-2; i++) {
    long start = (1L << i) - 1;
    long finish = (1L << (i+1)) - 2;
    for (long j = finish; j >= start; j--) {
      index_vec[2*j+2] = index_vec[j] + (index_vec[j+1] - index_vec[j])/2;
      index_vec[2*j+1] = index_vec[j];
    }
    index_vec[2*finish+3] = n;
  }
  
  for (long i = (1L << (levels-1)) - 1; i < vec_len; i++) {
    /* multiply primes index_vec[i]..index_vec[i+1]-1 into 
     * prod_vec[i]
     */
    _ntl_gone(&prod_vec[i]);
    for (long j = index_vec[i]; j < index_vec[i+1]; j++)
      _ntl_gsmul(prod_vec[i], q[j], &prod_vec[i]);
  }
  
  for (long i = (1L << (levels-1)) - 1; i < vec_len; i++) {
    for (long j = index_vec[i]; j < index_vec[i+1]; j++)
      _ntl_gsdiv(prod_vec[i], q[j], &coeff_vec[j]);
  }

  for (long i = (1L << (levels-1)) - 2; i >= 0; i--)
    _ntl_gmul(prod_vec[2*i+1], prod_vec[2*i+2], &prod_vec[i]);
  
  /*** new asymptotically fast code to compute inv_vec ***/
  
  _ntl_gone(&rem_vec[0]);
  for (long i = 0; i < (1L << (levels-1)) - 1; i++) {
    _ntl_gmod(rem_vec[i], prod_vec[2*i+1], &temps[0]);
    _ntl_gmul(temps[0], prod_vec[2*i+2], &temps[1]);
    _ntl_gmod(temps[1], prod_vec[2*i+1], &rem_vec[2*i+1]);
    
    _ntl_gmod(rem_vec[i], prod_vec[2*i+2], &temps[0]);
    _ntl_gmul(temps[0], prod_vec[2*i+1], &temps[1]);
    _ntl_gmod(temps[1], prod_vec[2*i+2], &rem_vec[2*i+2]);
  }

  for (long i = (1L << (levels-1)) - 1; i < vec_len; i++) {
    for (long j = index_vec[i]; j < index_vec[i+1]; j++) {
      long tt, tt1, tt2;
      _ntl_gsdiv(prod_vec[i], q[j], &temps[0]);
      tt = _ntl_gsmod(temps[0], q[j]);
      tt1 = _ntl_gsmod(rem_vec[i], q[j]);
      tt2 = MulMod(tt, tt1, q[j]);
      inv_vec[j] = sp_inv_mod(tt2, q[j]);
    }
  }
}

/*------------------------------------------------------------*/
/* Not applicable. Todo: make this cleaner                    */
/*------------------------------------------------------------*/
void ZZ_CRT_crt_fast::insert(long i, _ntl_gbigint m){
   LogicError("insert called improperly");
}

/*------------------------------------------------------------*/
/* CRT for ZZ_CRT_crt_fast                                    */
/*------------------------------------------------------------*/
void ZZ_CRT_crt_fast::crt(ZZ& x_ZZ, const Vec<long>& b_in){
  
  long *val_vec1 = val_vec._vec__rep;
  long i;
  
  for (i = 0; i < n; i++) {
    val_vec1[i] = MulMod(b_in[i], inv_vec[i], moduli[i]);
  }
  
  for (i = (1L << (levels-1)) - 1; i < vec_len; i++) {
    long j1 = index_vec[i];
    long j2 = index_vec[i+1];
    gadd_mul_many(&rem_vec[i], &coeff_vec[j1], &val_vec1[j1], j2-j1, SIZE(prod_vec[i]));
  }
  
  for (i = (1L << (levels-1)) - 2; i >= 0; i--) {
    _ntl_gmul(prod_vec[2*i+1], rem_vec[2*i+2], &temps[0]);
    _ntl_gmul(rem_vec[2*i+1], prod_vec[2*i+2], &temps[1]);
    _ntl_gadd(temps[0], temps[1], &rem_vec[i]);
  }
  
  /* temps[0] = rem_vec[0] mod prod_vec[0] (least absolute residue) */
  _ntl_gmod(rem_vec[0], prod_vec[0], &temps[0]);
  _ntl_gsub(temps[0], prod_vec[0], &temps[1]);
  _ntl_gnegate(&temps[1]);
  if (_ntl_gcompare(temps[0], temps[1]) > 0) {
    _ntl_gnegate(&temps[1]);
    _ntl_gcopy(temps[1], &temps[0]);
  }
  
  _ntl_gmod(temps[0], p.rep, &temps[1]);
  _ntl_gcopy(temps[1], &x_ZZ.rep);
}


/*------------------------------------------------------------*/
/* constructor of ZZ_CRT_crt_basic                            */
/* call to make at the end                                    */
/*------------------------------------------------------------*/
ZZ_CRT_crt_basic::ZZ_CRT_crt_basic(const Vec<long>& q){

  if (q.length() >= 600)
    return;

  moduli = q;
  n = moduli.length();

  p = to_ZZ(1);
  for (long i = 0; i < n; i++)
    p *= moduli[i];

  v.SetLength(n);
  for (long i = 0; i < n; i++)
    v[i] = 0;

  sbuf = SIZE(p.rep)+2;
  make();
}

/*------------------------------------------------------------*/
/* CRT for ZZ_CRT_crt_basic                                   */
/* premultiply at the beginning, reduce.eval at the end       */
/*------------------------------------------------------------*/
void ZZ_CRT_crt_basic::crt(ZZ& x_ZZ, const Vec<long>& b_in){
  Vec<long> b_vec;
  premultiply(b_vec, b_in);

  ZZ xx1;
  _ntl_gbigint* x = &xx1.rep;

   mp_limb_t *xx, *yy; 
   _ntl_gbigint x1;
   long i, sx;
   long sy;
   mp_limb_t carry;

   sx = sbuf;
   _ntl_gsetlength(x, sx);
   x1 = *x;
   xx = DATA(x1);

   for (i = 0; i < sx; i++)
      xx[i] = 0;

   for (i = 0; i < n; i++) {
      if (!v[i]) continue;

      yy = DATA(v[i]);
      sy = SIZE(v[i]); 

      if (!sy || !b_vec[i]) continue;

      carry = mpn_addmul_1(xx, yy, sy, b_vec[i]);
      yy = xx + sy;
      *yy += carry;

      if (*yy < carry) { /* unsigned comparison! */
         do {
            yy++;
            *yy += 1;
         } while (*yy == 0);
      }
   }

   while (sx > 0 && xx[sx-1] == 0) sx--;
   SIZE(x1) = sx;

   reduce_struct.eval(x_ZZ, xx1);
}

/*------------------------------------------------------------*/
/* inserts m into v[i]                                        */
/*------------------------------------------------------------*/
void ZZ_CRT_crt_basic::insert(long i, _ntl_gbigint m){
  _ntl_gcopy(m, &v[i]);
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* a class that does everything                               */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* helper routine to build rem attribute                      */
/*------------------------------------------------------------*/
void ZZ_CRT::build_rem_struct(const Vec<long>& moduli) {

  long n = moduli.length();
  
#if (defined(NTL_VIABLE_LL) && defined(NTL_TBL_REM))
  if (n <= NTL_MAX_REM_TBL){
    rem_struct = &rem_tbl;
    return;
  }
#endif
  if (n >= NTL_MIN_REM_MEDIUM && n < NTL_MIN_REM_FAST) {
    rem_struct = &rem_medium;
    return;
  }
  if (n >= NTL_MIN_REM_FAST) {
    rem_struct = &rem_fast;
    return;
  }
  else {
    rem_struct = &rem_basic;
    return;
  }
}

/*------------------------------------------------------------*/
/* helper routine to build crt attribute                      */
/*------------------------------------------------------------*/
void ZZ_CRT::build_crt_struct(const Vec<long>& moduli) {

  long n = moduli.length();

  if (n >= NTL_MAX_CRT_TBL){
    crt_struct = &crt_fast;
    return;
  }
#if (defined(NTL_VIABLE_LL)) && (defined(NTL_CRT_ALTCODE))
  {
    crt_struct = &crt_tbl;
    return;
  }
#endif
#if (defined(NTL_VIABLE_LL)) && (defined(NTL_CRT_ALTCODE_SMALL))
  if (n <= 16){
    crt_struct = &crt_tbl;
    return;
  }
  else{
    crt_struct = &crt_basic;
    return;
  }
#endif
  {
    crt_struct = &crt_basic;
    return;
  }
}


/*------------------------------------------------------------*/
/* constructor for ZZ_CRT                                     */
/*------------------------------------------------------------*/
ZZ_CRT::ZZ_CRT(const Vec<long>& moduli) : 
#if (defined(NTL_VIABLE_LL) && defined(NTL_TBL_REM))
  rem_tbl(moduli),
#endif
  rem_fast(moduli),
  rem_medium(moduli),
  rem_basic(moduli),
#if ((defined(NTL_VIABLE_LL)) && (defined(NTL_CRT_ALTCODE))) ||  ((defined(NTL_VIABLE_LL)) && (defined(NTL_CRT_ALTCODE_SMALL)))
  crt_tbl(moduli),
#endif
  crt_fast(moduli),
  crt_basic(moduli) {
  build_rem_struct(moduli);  
  build_crt_struct(moduli);  
  }

