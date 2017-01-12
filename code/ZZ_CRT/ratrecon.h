#ifndef RATRECON__H
#define RATRECON__H

#include <NTL/ZZ.h>
#include <NTL/vector.h>

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* rational reconstruction                                    */
/* Supports rational reconstruction of a Vec^i<ZZ> T          */
/* requires a witness = one Vec^i<long> and one prime p       */
/* such that T = a mod p                                      */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/


/*------------------------------------------------------------*/
/* computes a = b d mod m, for  Vec^i<ZZ>                     */
/* base case (i=0)                                            */
/*------------------------------------------------------------*/
static inline
void mulmod (ZZ& a, const ZZ& b, const ZZ& d, const ZZ& m){
  MulMod(a, b, d, m);
}

/*------------------------------------------------------------*/
/* computes a = b d mod m, for  Vec^i<ZZ>                     */
/* template for i > 0                                         */
/*------------------------------------------------------------*/
template <class T>
void mulmod(Vec<T>& a, const Vec<T>& b, const ZZ& d, const ZZ& m){
  long n = b.length();
  a.SetLength(n);
  for (long i = 0; i < n; i++)
    mulmod(a[i], b[i], d, m);
}


/*------------------------------------------------------------*/
/* a is a Vec^i<ZZ>                                           */
/* tests if all entries of a are less than bound in abs       */
/* base case (i=0)                                            */
/*------------------------------------------------------------*/
static inline
bool less_than(const ZZ& a, const ZZ& bound){
  return abs(a) < bound;
}

/*------------------------------------------------------------*/
/* a is a Vec^i<ZZ>                                           */
/* tests if all entries of a are less than bound in abs       */
/* template for i > 0                                         */
/*------------------------------------------------------------*/
template <class T>
bool less_than(const Vec<T>& a, const ZZ& bound){
  for (long i = 0; i < a.length(); i++)
    if (!less_than(a[i], bound))
      return false;
  return true;
}


/*------------------------------------------------------------*/
/* a is a Vec^i<ZZ> less than m                               */
/* subtracts m to all entries that are > m/2                  */
/* base case (i=0)                                            */
/*------------------------------------------------------------*/
static inline 
void center(ZZ& ac, const ZZ& a, const ZZ& m){
  if ((m-a) < a)
    ac = a-m;
  else
    ac = a;
}

/*------------------------------------------------------------*/
/* a is a Vec^i<ZZ> less than m                               */
/* subtracts m to all entries that are > m/2                  */
/* template for i > 0                                         */
/*------------------------------------------------------------*/
template <class T>
void center(Vec<T>& ac, const Vec<T>& a, const ZZ& m){
  ac.SetLength(a.length());
  for (long i = 0; i < a.length(); i++)
    center(ac[i], a[i], m);
}


/*------------------------------------------------------------*/
/* computes a = b d, where b is a Vec^i<ZZ>                   */
/* template for i > 0 (no need to write i=0; built-in)        */
/*------------------------------------------------------------*/
template <class T>
void mul(Vec<T>& a, const Vec<T>& b, const ZZ& d){
  long n = b.length();
  a.SetLength(n);
  for (long i = 0; i < n; i++)
    mul(a[i], b[i], d);
}


/*------------------------------------------------------------*/
/* a is a Vec^i<ZZ>, witness is a Vec^i<long>                 */
/* tests whether a/denom = witness mod prime                  */
/* base case (i=0)                                            */
/*------------------------------------------------------------*/
static inline
bool check_mod(const ZZ& a, const long denom, const long& witness, const long prime){
  return (a % prime) == MulMod(denom, witness, prime);
}

/*------------------------------------------------------------*/
/* a is a Vec^i<ZZ>, witness is a Vec^i<long>                 */
/* tests whether a/denom = witness mod prime                  */
/* template for i > 0                                         */
/*------------------------------------------------------------*/
template <class T, class U>
bool check_mod(const Vec<T>& a, const long denom, const Vec<U>& witness, const long prime){
  long n = a.length();
  if (witness.length() != n)
    return false;
  for (long i = 0; i < n; i++){
    bool b = check_mod(a[i], denom, witness[i], prime);
    if (b == false)
      return false;
  }
  return true;
}


/*-----------------------------------------------------------*/
/* rational reconstruction                                   */
/*   x is a Vec^i<ZZ> with entries in [0..m-1]               */
/*   witness is a Vec^i<long> of the same size               */
/*   half = m^(1/2), most=m^(3/4)                            */
/* checks whether                                            */
/*    x reconstructs to a/b and                              */
/*    a/b = witness*multiplier mod prime                     */
/* base case (i=0); uses NTL's rational reconstruction       */
/*    most=m^(3/4) is not used here                          */
/*-----------------------------------------------------------*/
static inline
long ReconstructRational_and_check(ZZ& a, ZZ& b, const ZZ& x, const ZZ& m, const ZZ& half, const ZZ& most, 
				   const long& witness, const long prime,  long multiplier){
  long s = ReconstructRational(a, b, x, m, half, half);
  if (s == 0)
    return 0;
  if ((a % prime) == MulMod(MulMod(b % prime, witness, prime), multiplier, prime))
    return 1;
  else 
    return 0;
}

/*-----------------------------------------------------------*/
/* rational reconstruction                                   */
/*   x is a Vec^i<ZZ> with entries in [0..m-1]               */
/*   witness is a Vec^i<long> of the same size               */
/*   half = m^(1/2), most=m^(3/4)                            */
/* checks whether                                            */
/*    x reconstructs to a/b and                              */
/*    a/b = witness*multiplier mod prime                     */
/* template (i > 0) uses previous denominators:              */
/*    we pre-multiply u[i] by the denoms found before        */
/*    if the result is "small" (less than m^(3/4) in abs)    */
/*    we keep it; otherwise we do a rational reconstruction  */
/*-----------------------------------------------------------*/
template <class T, class U>
long ReconstructRational_and_check(Vec<T>& a, ZZ& b, const Vec<T>& u, const ZZ& m, const ZZ& half_prec, const ZZ& most_prec,
				   const Vec<U>& witness, const long prime, long multiplier){
  long n = u.length();
  a.SetLength(n);
  
  Vec<ZZ> dens, prod, inv_prod;
  prod.SetLength(n+1); // product of all denominators found so far
  prod[n] = to_ZZ(1);
  dens.SetLength(n); // the new denominator found at step i
  
  Vec<long> prod_mod;
  prod_mod.SetLength(n+1); // product of all denominators found so far, times the initial multiplier, all taken modulo the new prime
  prod_mod[n] = multiplier;
  
  Vec<T> tmp_a1, tmp_a1_center, tmp_a2;
  tmp_a1.SetLength(n);  // tmp_a1[i] = u[i] * (all denominators seen so far)
  tmp_a1_center.SetLength(n); // the same, after re-centering the entries
  tmp_a2.SetLength(n);  // tmp_a2[i] = numerators of rat-recon of tmp_a1[i]

  for (long i = n-1; i >= 0; i--){
    mulmod(tmp_a1[i], u[i], prod[i+1], m);
    center(tmp_a1_center[i], tmp_a1[i], m);

    if (less_than(tmp_a1_center[i], most_prec)){ // if all entries of tmp_a1_centered are < m^(3/4) in absolute value
                                                 // looks like we have found an integer
      tmp_a2[i] = tmp_a1_center[i];              // -> no need to rat-recon anything here
      bool boo = check_mod(tmp_a2[i], prod_mod[i+1], witness[i], prime); // check it with the witness
      if (boo == false)
	return 0;
      dens[i] = to_ZZ(1);           // if OK, no new denominator
      prod[i] = prod[i+1];          // product stays the same
      prod_mod[i] = prod_mod[i+1];  // product modulo new prime too
    }
    else{
      long s = ReconstructRational_and_check(tmp_a2[i], dens[i], tmp_a1[i], m, half_prec, most_prec, witness[i], prime, prod_mod[i+1]);
      // checks whether tmp_a2[i] reconstructs
      // in the recursive call, we should multiply by (prod[i+1] mod prime), since we did the same to u[i]
      // this is done by using prod_mod[i+1] as a multiplier
      if (s == 0)
	return 0;
      prod[i] = dens[i]*prod[i+1];
      prod_mod[i] = MulMod(dens[i] % prime, prod_mod[i+1], prime);
    }
  }

  // now use the last denominator (=prod[0]) as a common denominator for every entry
  ZZ fac = to_ZZ(1);
  a[0] = tmp_a2[0];
  for (long i = 1; i < n; i++){
    fac = fac*dens[i-1];
    mul(a[i], tmp_a2[i], fac);
  }
  b = prod[0];
  return 1;
}


/*-----------------------------------------------------------*/
/* rational reconstruction                                   */
/*   x is a Vec^i<ZZ> with entries in [0..m-1]               */
/*   witness is a Vec^i<long> of the same size               */
/* checks whether                                            */
/*    x reconstructs to a/b and                              */
/*    a/b = witness            mod prime                     */
/*-----------------------------------------------------------*/
template <class T, class U>
  long ReconstructRational(Vec<T>& a, ZZ& b, const Vec<T>& u, const ZZ& m, const Vec<U>& witness, const long prime){
  ZZ half = SqrRoot(m/2);
  ZZ most = half*SqrRoot(half);

  return ReconstructRational_and_check(a, b, u, m, half, most, witness, prime, 1L);
}

#endif
