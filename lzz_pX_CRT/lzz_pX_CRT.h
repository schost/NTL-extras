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
  /* constructor for geometric progressione                     */
  /* we interpolate at a^(2(i+i0)), i=0..N-1                    */
  /* assume all polynomials we evaluate have degree < eval_N    */
  /*------------------------------------------------------------*/
  zz_pX_Multipoint_Geometric(const zz_p& a, long eval_N, long N, long i0);

  /*------------------------------------------------------------*/
  /* basic operations                                           */
  /*------------------------------------------------------------*/
  void evaluate_one(Vec<zz_p>& val, const zz_pX& f) const;
  void evaluate(Vec<zz_p>& val, const zz_pX& f) const;
  void evaluate(Vec<Vec<zz_p>>& val, const Vec<zz_pX>& f) const;
  void interpolate(zz_pX& f, const Vec<zz_p>& val) const;

  
  ~zz_pX_Multipoint_Geometric(){};

  // private:
  long eval_n;
  long i0;
  Vec<zz_p> eval_inverse_powers_square_q, eval_inverse_powers_square_q_shifted, inverse_powers_square_q, inverse_powers_square_q_shifted, inverse_derivative;
  zz_pX M, revM, eval_S, S;
  zz_p q;
  fftRep eval_S_FFT, S_FFT, revM_FFT;
  Mat<zz_p> shift_matrix;
};



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

/*------------------------------------------------------------*/
/* multiplicative order of a                                  */
/* -1 if a = 0                                                */
/*------------------------------------------------------------*/
long order(const zz_p& a);

/*------------------------------------------------------------*/
/* finds an element of order at least ord                     */
/* assumes it exists, does not verify                         */
/*------------------------------------------------------------*/
void element_of_order(zz_p& a, long ord);

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* Chinese Remaindernig over zz_p                             */
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
