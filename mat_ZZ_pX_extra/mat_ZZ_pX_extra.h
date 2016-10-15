#ifndef MAT_ZZ_PX_EXTRA__H
#define MAT_ZZ_PX_EXTRA__H

#include <NTL/matrix.h>
#include <NTL/ZZ_pX.h>

NTL_CLIENT

/*------------------------------------------------------------*/
/* magma output, without assign, using variable var           */
/*------------------------------------------------------------*/
void magma_output(const Mat<ZZ_pX>& a, const string & var);

/*------------------------------------------------------------*/
/* magma assignment, using variable var, to name              */
/*------------------------------------------------------------*/
void magma_assign(const Mat<ZZ_pX>& a, const string & var, const string & name);

/*------------------------------------------------------------*/
/* random matrix of a given degree                            */
/*------------------------------------------------------------*/
void random_mat_ZZ_pX(Mat<ZZ_pX>& a, long n, long m, long d);

/*------------------------------------------------------------*/
/* maximum degree of the entries of a                         */
/*------------------------------------------------------------*/
long deg(const Mat<ZZ_pX> & a);

/*------------------------------------------------------------*/
/* maximum alloc of the entries of a                          */
/*------------------------------------------------------------*/
long alloc(const Mat<ZZ_pX>& a);

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/*  direct multiplication                                     */
/*------------------------------------------------------------*/
void mul_direct(Mat<ZZ_pX> & c, const Mat<ZZ_pX> & a, const Mat<ZZ_pX> & b);

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/*  Wakman's algorithm                                        */
/*------------------------------------------------------------*/
void mul_waksman(Mat<ZZ_pX> & c, const Mat<ZZ_pX> & a, const Mat<ZZ_pX> & b);

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/*  reduction to a multiplication modulo a small FFT prime    */
/*------------------------------------------------------------*/
void mul_CRT_CTFT(Mat<ZZ_pX>& C, const Mat<ZZ_pX>& A, const Mat<ZZ_pX>& B);

#endif
