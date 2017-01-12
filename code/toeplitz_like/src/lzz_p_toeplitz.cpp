#include <NTL/lzz_pX.h>
#include <NTL/mat_lzz_p.h>

#include "lzz_pX_middle_product.h"
#include "lzz_p_toeplitz.h"

NTL_CLIENT

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/* Toeplitz matrices                                  */
/* stored as                                          */
/*       a2 a3 a4 a5                                  */
/*       a1 a2 a3 a4                                  */
/*       a0 a1 a2 a3                                  */
/*----------------------------------------------------*/
/*----------------------------------------------------*/

/*----------------------------------------------------*/
/* default constructor                                */
/*----------------------------------------------------*/
lzz_p_toeplitz::lzz_p_toeplitz(){
  data.SetLength(0);
  n = m = 0;
}

/*----------------------------------------------------*/
/* input vector is as showed above                    */
/*----------------------------------------------------*/
lzz_p_toeplitz::lzz_p_toeplitz(const Vec<zz_p>& input, long rows, long cols){
  n = rows;
  m = cols;
  data = input;
  data_rev.SetLength(n+m-1);
  zz_pX data_X;
  data_X.rep.SetLength(n+m-1);
  
  for (long i = 0;i < n+m-1; i++){
    data_rev[i] = input[n+m-2-i];
    data_X.rep[i] = data_rev[i];
  }
  data_X.normalize();
  
  long K = NextPowerOfTwo(n+m-1);
  fft_data = fftRep(INIT_SIZE, K);
  TofftRep(fft_data, data_X, K);

  // TODO: do not init if not FFT prime
  c = zz_pX_Multipoint_CTFT(n+m-1);
  c.interpolate_t(TFT_data, input);
}

/*---------------------------------------------------*/
/* dimensions                                        */
/*---------------------------------------------------*/
long lzz_p_toeplitz::NumRows() const {
  return n;
}

long lzz_p_toeplitz::NumCols() const {
  return m;
}

/*---------------------------------------------------*/
/* data access                                       */
/*---------------------------------------------------*/
const zz_p& lzz_p_toeplitz::operator ()(long i, long j) const {
  return data[n-1-i+j];
}

/*----------------------------------------------------*/
/* turns M into a dense matrix                        */
/*----------------------------------------------------*/
void lzz_p_toeplitz::to_dense(Mat<zz_p>& Mdense) const { 
  long n = NumRows();
  long m = NumCols();
  Mdense.SetDims(n, m);
  
  for (long i = 0; i < n; i++)
    for (long j = 0; j < m; j++)
      Mdense[i][j] = (*this)(i,j);
}

/*----------------------------------------------------*/
/* right multiplication                               */
/*----------------------------------------------------*/
void lzz_p_toeplitz::mul_right(Vec<zz_p>& res, const Vec<zz_p>& input) const {

  long nM = NumRows();
  long mM = NumCols();
  res.SetLength(nM);

  if (nM == 0 || mM == 0){
    for (long i = 0; i < nM; i++)
      res[i] = 0;
    return;
  }
  if (min(nM, mM) <= NTL_zz_pX_MUL_CROSSOVER / 2){
    Vec<zz_p> input_rev;
    input_rev.SetLength(mM);
    zz_p *cf = input_rev.elts();
    for (long i = 0; i < mM; i++)
      cf[i] = input[mM-1-i];

    long sp = Kar_stk_size(max(nM, mM));
    zz_p *stk = new zz_p[sp];
    tKarMul_aux(res._vec__rep.rep, nM, input_rev._vec__rep.rep, mM, data_rev._vec__rep.rep, nM+mM-1, stk);
    delete[] stk;
  }
  else {
    long K = NextPowerOfTwo(nM+mM-1);
    fftRep fft_input = fftRep(INIT_SIZE, K);

    zz_pX input_X;
    input_X.rep.SetLength(mM);
    zz_p *cf = input_X.rep.elts();
    for (long i = 0; i < mM; i++)
      cf[i] = input[i];
    input_X.normalize();

    TofftRep(fft_input, input_X, K);
    mul(fft_input, fft_input, fft_data);
    FromfftRep(res._vec__rep.rep, fft_input, mM-1, nM+mM-2);
  }
}

/*----------------------------------------------------*/
/* right multiplication                               */
/*----------------------------------------------------*/
void lzz_p_toeplitz::mul_right_CTFT(Vec<zz_p>& res, const Vec<zz_p>& input) const {

  long m = NumRows();
  long n = NumCols();
  res.SetLength(m);

  zz_pX input_f;
  input_f.rep = input;
  input_f.normalize();

  Vec<zz_p> TFT_input, res_rev;
  c.evaluate(TFT_input, input_f);
  for (long i = 0; i < m+n-1; i++)
    TFT_input[i] *= TFT_data[i];
  
  c.evaluate_t(res_rev, TFT_input);
  res.SetLength(m);
  for (long i = 0; i < m; i++)
    res[i] = res_rev[m-1-i];
}

/*----------------------------------------------------*/
/* right multiplication                               */
/*----------------------------------------------------*/
void lzz_p_toeplitz::mul_right(Mat<zz_p>& res, const Mat<zz_p>& input) const {
  long m = NumRows();
  long n = NumCols();

  Vec<zz_p> vec_in, vec_out;
  vec_in.SetLength(n);
  long a = input.NumCols();
  res.SetDims(m, a);
  for (long i = 0; i < a; i++){
    zz_p *elts = vec_in.elts();
    for (long j = 0; j < n; j++)
      elts[j] = input[j][i];
    mul_right(vec_out, vec_in);
    elts = vec_out.elts();
    for (long j = 0; j < m; j++)
      res[j][i] = elts[j];
  }
}

/*----------------------------------------------------*/
/* left multiplication                                */
/*----------------------------------------------------*/
void lzz_p_toeplitz::mul_left(Vec<zz_p>& res, const Vec<zz_p>& input) const {

  long nM = NumRows();
  long mM = NumCols();

  res.SetLength(mM);


  if (nM == 0 || mM == 0){
    for (long i = 0; i < mM; i++)
      res[i] = 0;
    return;
  }
  if (min(nM, mM) <= NTL_zz_pX_MUL_CROSSOVER / 2){
    Vec<zz_p> output_rev;
    output_rev.SetLength(mM);
    res.SetLength(mM);

    long sp = Kar_stk_size(max(nM, mM));
    zz_p *stk = new zz_p[sp];
    tKarMul_aux(output_rev._vec__rep.rep, mM, input._vec__rep.rep, nM, data_rev._vec__rep.rep, nM+mM-1, stk);
    delete[] stk;

    zz_p *cf = output_rev.elts();
    zz_p *cf_res = res.elts();
    for (long i = 0; i < mM; i++)
      cf_res[i] = cf[mM-1-i];
  }
  else{
    long K = NextPowerOfTwo(nM+mM-1);
    fftRep fft_input = fftRep(INIT_SIZE, K);

    zz_pX input_X;
    input_X.rep.SetLength(nM);
    zz_p *cf = input_X.rep.elts();
    for (long i = 0; i < nM; i++)
      cf[i] = input[nM-1-i];
    input_X.normalize();

    Vec<zz_p> output_rev;
    output_rev.SetLength(mM);

    TofftRep(fft_input, input_X, K);
    mul(fft_input, fft_input, fft_data);
    FromfftRep(output_rev._vec__rep.rep, fft_input, nM-1, nM+mM-2);

    res.SetLength(mM);
    cf = output_rev.elts();
    zz_p *cf_res = res.elts();
    for (long i = 0; i < mM; i++)
      cf_res[i] = cf[mM-1-i];
  }
}

/*----------------------------------------------------*/
/* left multiplication                                */
/*----------------------------------------------------*/
void lzz_p_toeplitz::mul_left(Mat<zz_p>& res, const Mat<zz_p>& input) const {
  long m = NumRows();
  long n = NumCols();

  Vec<zz_p> vec_in, vec_out;
  vec_in.SetLength(m);
  long a = input.NumCols();
  res.SetDims(n, a);
  for (long i = 0; i < a; i++){
    zz_p *elts = vec_in.elts();
    for (long j = 0; j < m; j++)
      elts[j] = input[j][i];
    mul_left(vec_out, vec_in);
    elts = vec_out.elts();
    for (long j = 0; j < n; j++)
      res[j][i] = elts[j];
  }
}
