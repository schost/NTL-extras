#ifndef MAGMA_OUTPUT__H
#define MAGMA_OUTPUT__H

#include <NTL/ZZX.h>
#include <NTL/lzz_pX.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* basic functions for output in magma format                 */
/* handles mainly Vec<zz_p>, zz_pX, ZZX                       */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* prints a matrix                                            */
/*------------------------------------------------------------*/
void magma_output(const Mat<zz_p> & v);

/*------------------------------------------------------------*/
/* assigns a matrix to variable "name"                        */
/*------------------------------------------------------------*/
void magma_assign(const Mat<zz_p> & v, const string & name);

/*------------------------------------------------------------*/
/* prints a vector                                            */
/*------------------------------------------------------------*/
void magma_output(const Vec<ZZ> & v);

/*------------------------------------------------------------*/
/* assigns a vector to variable "name"                        */
/*------------------------------------------------------------*/
void magma_assign(const Vec<ZZ> & v, const string & name);

/*------------------------------------------------------------*/
/* initializes UZ<XX>=Z[XX]                                   */
/*------------------------------------------------------------*/
void magma_init_ZZX();

/*------------------------------------------------------------*/
/* prints a poly with indeterminate var                       */
/*------------------------------------------------------------*/
void magma_output(const ZZX & v, const string & var);

/*------------------------------------------------------------*/
/* prints a poly with indeterminate x, as a cast into U       */
/*------------------------------------------------------------*/
void magma_output(const ZZX & v);

/*------------------------------------------------------------*/
/* assign a poly with indeterminate var to variable "name"    */
/*------------------------------------------------------------*/
void magma_assign(const ZZX & v, const string & var, const string & name);

/*------------------------------------------------------------*/
/* assign a poly with indeterminate x to variable "name"      */
/*------------------------------------------------------------*/
void magma_assign(const ZZX & v, const string & name);

/*------------------------------------------------------------*/
/* prints a vector                                            */
/*------------------------------------------------------------*/
void magma_output(const Vec<long> & v);

/*------------------------------------------------------------*/
/* assigns a vector to variable "name"                        */
/*------------------------------------------------------------*/
void magma_assign(const Vec<long> & v, const string & name);

/*------------------------------------------------------------*/
/* initializes GF(p)                                          */
/*------------------------------------------------------------*/
void magma_init();

/*------------------------------------------------------------*/
/* initializes U<x>=GF(p)[x]                                  */
/*------------------------------------------------------------*/
void magma_init_X();

/*------------------------------------------------------------*/
/* prints a vector                                            */
/*------------------------------------------------------------*/
void magma_output(const Vec<zz_p> & v);

/*------------------------------------------------------------*/
/* assigns a vector to variable "name"                        */
/*------------------------------------------------------------*/
void magma_assign(const Vec<zz_p> & v, const string & name);

/*------------------------------------------------------------*/
/* prints a poly with indeterminate var                       */
/*------------------------------------------------------------*/
void magma_output(const zz_pX & v, const string & var);

/*------------------------------------------------------------*/
/* prints a poly with indeterminate x, as a cast into U       */
/*------------------------------------------------------------*/
void magma_output(const zz_pX & v);

/*------------------------------------------------------------*/
/* assign a poly with indeterminate var to variable "name"    */
/*------------------------------------------------------------*/
void magma_assign(const zz_pX & v, const string & var, const string & name);

/*------------------------------------------------------------*/
/* assign a poly with indeterminate x to variable "name"      */
/*------------------------------------------------------------*/
void magma_assign(const zz_pX & v, const string & name);


/*------------------------------------------------------------*/
/* prints a vector of polys with indeterminate var            */
/*------------------------------------------------------------*/
void magma_output(const Vec<zz_pX> & v, const string & var);

/*------------------------------------------------------------*/
/* prints a vector of polys with indeterminate x              */
/*------------------------------------------------------------*/
void magma_output(const Vec<zz_pX> & v);

/*------------------------------------------------------------*/
/* assign a vector of polys with indet var to variable "name" */
/*------------------------------------------------------------*/
void magma_assign(const Vec<zz_pX> & v, const string & var, const string & name);

/*------------------------------------------------------------*/
/* assign a vector of polys with indet x to variable "name"   */
/*------------------------------------------------------------*/
void magma_assign(const Vec<zz_pX> & v, const string & name);


#endif
