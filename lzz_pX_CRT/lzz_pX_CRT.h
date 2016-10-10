#ifndef LZZ_PX_CRT__H
#define LZZ_PX_CRT__H

#include <NTL/lzz_pX.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* classes for CRT and multipoint evaluation over zz_p        */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/


/*------------------------------------------------------------*/
/* builds the subproduct tree naively (no Huffman)            */
/*------------------------------------------------------------*/
void build_subproduct_tree(Vec<Vec<zz_pX>> & tree, const Vec<zz_pX> & q);


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* multipoint evaluation over zz_p                            */
/* abstract class                                             */
/* TODO: thresholds                                           */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class zz_pX_Multipoint{
 public:

  /*------------------------------------------------------------*/
  /* basic operations                                           */
  /*------------------------------------------------------------*/
  virtual void evaluate(Vec<zz_p>& val, const zz_pX& f) const = 0;
  virtual void evaluate(Vec<Vec<zz_p>>& val, const Vec<zz_pX>& f) const = 0;
  virtual void interpolate(zz_pX& f, const Vec<zz_p>& val) const = 0;
 
  /*------------------------------------------------------------*/
  /* number of points                                           */
  /*------------------------------------------------------------*/
  long length() const{
    return n;
  }

  /*------------------------------------------------------------*/
  /* points                                                     */
  /*------------------------------------------------------------*/
  void point(zz_p& pt, long i) const{
    pt = pts[i];
  }

  // TODO: why isn't this private?
  Vec<zz_p> pts;
  long n;
};


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* multipoint evaluation over zz_p                            */
/* uses subproduct tree.                                      */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class zz_pX_Multipoint_General : public zz_pX_Multipoint{
 public:
  
  zz_pX_Multipoint_General(){}

  /*------------------------------------------------------------*/
  /* constructor from points                                    */
  /*------------------------------------------------------------*/
  zz_pX_Multipoint_General(const Vec<zz_p>& q){
    pts = q;
    n = q.length();
    Vec<zz_pX> qpol;
    qpol.SetLength(n);
    for (long i = 0; i < n; i++){
      SetCoeff(qpol[i], 0, -q[i]);
      SetCoeff(qpol[i], 1, 1);
    }
    build_subproduct_tree(tree, qpol);
    
    evaluate(cofactors, diff(tree[tree.length()-1][0]));
    for (long i = 0; i < n; i++)
      cofactors[i] = 1/cofactors[i];
  }

  /*------------------------------------------------------------*/
  /* basic operations                                           */
  /*------------------------------------------------------------*/
  void evaluate(Vec<zz_p>& val, const zz_pX& f) const;
  void evaluate(Vec<Vec<zz_p>>& val, const Vec<zz_pX>& f) const;
  void interpolate(zz_pX& f, const Vec<zz_p>& val) const;
  
  ~zz_pX_Multipoint_General(){};

 private:
  Vec<zz_p> cofactors;
  Vec<Vec<zz_pX> > tree;
};


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* multipoint evaluation over zz_p                            */
/* points in geometric progression                            */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class zz_pX_Multipoint_Geometric : public zz_pX_Multipoint{
 public:

  zz_pX_Multipoint_Geometric(){}

  /*------------------------------------------------------------*/
  /* constructor for geometric progressions                     */
  /* we interpolate at a^(2(i+i0)), i=0..N-1                    */
  /*------------------------------------------------------------*/
  zz_pX_Multipoint_Geometric(const zz_p& a, long N, long i0);

  /*------------------------------------------------------------*/
  /* constructor for geometric progressions                     */
  /* we interpolate at c.a^(2i), i=0..N-1                       */
  /*------------------------------------------------------------*/
  zz_pX_Multipoint_Geometric(const zz_p& a, long N, const zz_p& c);

  /*------------------------------------------------------------*/
  /* constructor for geometric progressione                     */
  /* we interpolate at a^(2(i+i0)), i=0..N-1                    */
  /* assume all polynomials we evaluate have degree < eval_N    */
  /*------------------------------------------------------------*/
  zz_pX_Multipoint_Geometric(const zz_p& a, long eval_N, long N, long i0);

  /*------------------------------------------------------------*/
  /* basic operations                                           */
  /*------------------------------------------------------------*/
  void evaluate(Vec<zz_p>& val, const zz_pX& f) const;
  void evaluate(Vec<Vec<zz_p>>& val, const Vec<zz_pX>& f) const;
  void interpolate(zz_pX& f, const Vec<zz_p>& val) const;

  
  ~zz_pX_Multipoint_Geometric(){};

  // TODO: any reason they are not private?
  // private:  

  long eval_n;
  long i0;
  Vec<zz_p> eval_inverse_powers_square_q, eval_inverse_powers_square_q_shifted, inverse_powers_square_q, inverse_powers_square_q_shifted, inverse_derivative, powers_c, inv_powers_c;
  zz_pX M, revM, eval_S, S;
  zz_p q, c;
  fftRep eval_S_FFT, S_FFT, revM_FFT;
  Mat<zz_p> shift_matrix;

 private:
  void evaluate_one(Vec<zz_p>& val, const zz_pX& f) const;
};


/*------------------------------------------------------------*/
/* a naive conversion to a dense matrix                       */
/* maybe promote to all multipoint matrices?                  */
/*------------------------------------------------------------*/
void to_dense(Mat<zz_p>& M, const zz_pX_Multipoint_Geometric& X);


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* multipoint evaluation over zz_p                            */
/* FFT points                                                 */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class zz_pX_Multipoint_FFT : public zz_pX_Multipoint{
 public:

  zz_pX_Multipoint_FFT(){}

  zz_pX_Multipoint_FFT(long n){
    k = NextPowerOfTwo(n);
    if (n != (1L << k)){
      cerr << "Wrong length for zz_pX_Multipoint_FFT\n";
      exit(-1);
    }

    this->n = n;
    if (zz_pInfo->p_info == NULL){
      cerr << "Attempt to init a zz_pX_Multipoint_FFT without zz_p::FFTInit\n";
      exit(-1);
    }

    pts.SetLength(n);
    zz_pX X;
    SetCoeff(X, 1, to_zz_p(1));
    evaluate(pts, X);
  }

  /*------------------------------------------------------------*/
  /* basic operations                                           */
  /*------------------------------------------------------------*/
  void evaluate(Vec<zz_p>& val, const zz_pX& f) const;
  void evaluate(Vec<Vec<zz_p>>& val, const Vec<zz_pX>& f) const;
  void interpolate(zz_pX& f, const Vec<zz_p>& val) const;
  
  ~zz_pX_Multipoint_FFT(){};

 private:
  long k;
};






/*---------------------------------------------------------*/
/*---------------------------------------------------------*/
/* a class that contains precomputed information           */
/* for negacyclic convolution                              */
/* z, z_precomp are multipliers for forward transform      */
/* invz, invz_precomp are for the inverse transform        */
/*---------------------------------------------------------*/
/*---------------------------------------------------------*/

class zz_pX_Multipoint_CTFT : public zz_pX_Multipoint {
 public:

  /*------------------------------------------------------------*/
  /* basic operations                                           */
  /*------------------------------------------------------------*/
  void evaluate(Vec<zz_p>& val, const zz_pX& f) const;
  void evaluate(Vec<Vec<zz_p>>& val, const Vec<zz_pX>& f) const {}
  void interpolate(zz_pX& f, const Vec<zz_p>& val) const;
  
  ~zz_pX_Multipoint_CTFT(){};
 
  /*--------------------------------------------------------------*/
  /* initializes the k-th row of pre-multipliers for negacyclic   */
  /* convolution (== in size 2^k)                                 */
  /*--------------------------------------------------------------*/
  void init_multipliers(const long k);
  
  /*-----------------------------------------------------------*/
  /* applies the multi-reduction map to f (mod p), in place    */
  /* input is in [0,p), output is in [0,2p)                    */
  /* x must have size at least 2n, and upper half must be zero */
  /*-----------------------------------------------------------*/
  void reduce(long* f) const;
 
  /*-----------------------------------------------------------*/
  /* applies the CRT map to x (mod p), in place                */
  /* x has length n                                            */
  /* input is in [0,p), output is in [0,p)                     */
  /* tmp is a temporary workspace, of length at least 2n       */
  /*-----------------------------------------------------------*/
  void CRT(long* x, long *tmp) const;


  zz_pX_Multipoint_CTFT(){}
  zz_pX_Multipoint_CTFT(long index_in){
    MaxK = -1;
    n = index_in;
    init_multipliers(NextPowerOfTwo(n));
    init_points();
  }


 private:
  long MaxK;
  Vec< Vec<long> > z;
  Vec< Vec<mulmod_precon_t> > z_precomp;
  Vec< Vec<long> > invz;
  Vec< Vec<mulmod_precon_t> > invz_precomp;

  /*-----------------------------------------------------------*/
  /* sets all entries in the vector pts                        */ 
  /*-----------------------------------------------------------*/
  void init_points();

  /*-----------------------------------------------------------*/
  /* evaluates a chunk of length 2^k                           */
  /* input length = 2^k                                        */
  /* input is word-size, output is in [0,p)                    */
  /* wk is a workspace of length 2^k                           */
  /*-----------------------------------------------------------*/
  void evaluate_chunk(zz_p * val, const long * coeffs, const long k, long * wk) const;
 
 /*------------------------------------------------------------*/
  /* interpolates a chunk of length 2^k                        */
  /* input length = 2^k                                        */
  /* wk is a workspace of length 2^k                           */
  /*-----------------------------------------------------------*/
  void interpolate_chunk(long * val, const zz_p * coeffs, const long k, long * wk) const;
 
  /*-----------------------------------------------------------*/
  /* applies postmultiplier after negacyclic convolution       */
  /* input length = 2^k                                        */
  /* input is word-size, output is in [0,p)                    */
  /* incorporates the division by 2^k for FFT^(-1)             */
  /*-----------------------------------------------------------*/
  void rescale_inverse(long *y, const long* x, const long k);
};

/*--------------------------------------------------------------*/
/* returns the exponents ki such that                           */
/*    n = sum_i 2^k_i                                           */
/*--------------------------------------------------------------*/
void CTFT_exponents(Vec<long>& expo, const long nn);




/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* Chinese Remaindering over zz_p                             */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class zz_pX_CRT{
 public:
  
  /*------------------------------------------------------------*/
  /* constructor from moduli                                    */
  /*------------------------------------------------------------*/
  zz_pX_CRT(const Vec<zz_pX>& q){
    n = q.length();
    build_subproduct_tree(tree, q);
    
    multimod(cofactors, diff(tree[tree.length()-1][0]));
    for (long i = 0; i < n; i++)
      cofactors[i] = MulMod(InvMod(cofactors[i], tree[0][i]), diff(tree[0][i]), tree[0][i]);
  }
	
  ~zz_pX_CRT(){};
    
  /*------------------------------------------------------------*/
  /* basic operations                                           */
  /*------------------------------------------------------------*/
  void multimod(Vec<zz_pX>& val, const zz_pX& f) const;
  void combine(zz_pX& f, const Vec<zz_pX>& val) const;
 
  /*------------------------------------------------------------*/
  /* master polynomial                                          */
  /*------------------------------------------------------------*/
  inline zz_pX master() const{
    return tree[tree.length()-1][0];
  }

  /*------------------------------------------------------------*/
  /* access to moduli                                           */
  /*------------------------------------------------------------*/
  inline zz_pX moduli(long i) const{
    return tree[0][i];
  }

  /*------------------------------------------------------------*/
  /* number of moduli                                           */
  /*------------------------------------------------------------*/
  long length() const{
    return n;
  }

 private:
  long n;
  Vec<zz_pX> cofactors;
  Vec<Vec<zz_pX>> tree;

};

#endif
