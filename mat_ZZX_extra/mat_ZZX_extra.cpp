#ifndef MAT_ZZX_EXTRA__H
#define MAT_ZZX_EXTRA__H

#include <NTL/matrix.h>
#include <NTL/lzz_pX.h>

NTL_CLIENT

/* /\*------------------------------------------------------------*\/ */
/* /\* magma output, without assign, using variable var           *\/ */
/* /\*------------------------------------------------------------*\/ */
/* void magma_output(const Mat<ZZX>& a, const string & var); */

/* /\*------------------------------------------------------------*\/ */
/* /\* magma assignment, using variable var, to name              *\/ */
/* /\*------------------------------------------------------------*\/ */
/* void magma_assign(const Mat<ZZX>& a, const string & var, const string & name); */

/* /\*------------------------------------------------------------*\/ */
/* /\* random matrix of a given degree                            *\/ */
/* /\*------------------------------------------------------------*\/ */
/* void random_mat_ZZX(Mat<ZZX>& a, long n, long m, long d); */

/* /\*------------------------------------------------------------*\/ */
/* /\* maximum degree of the entries of a                         *\/ */
/* /\*------------------------------------------------------------*\/ */
/* long deg(const Mat<ZZX> & a); */

/* /\\*------------------------------------------------------------*\/ */
/* /\* c = a*b                                                    *\/ */
/* /\*------------------------------------------------------------*\/ */
/* void multiply_waksman(Mat<ZZX> & c, const Mat<ZZX> & a, const Mat<ZZX> & b); */
/* void multiply_naive(Mat<ZZX> & c, const Mat<ZZX> & a, const Mat<ZZX> & b); */
/* void multiply_evaluate(Mat<ZZX> & c, const Mat<ZZX> & a, const Mat<ZZX> & b); */



#endif
