#ifndef SAGE_OUTPUT__H
#define SAGE_OUTPUT__H

#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* basic functions for output in sage format                  */
/* handles Vec<zz_p> and zz_pX                                */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/


/*------------------------------------------------------------*/
/* initializes GF(p)                                          */
/*------------------------------------------------------------*/
void sage_init();

/*------------------------------------------------------------*/
/* initializes U.<x>=GF(p)[x]                                 */
/*------------------------------------------------------------*/
void sage_init_X();

/*------------------------------------------------------------*/
/* prints a vector                                            */
/*------------------------------------------------------------*/
void sage_output(const Vec<zz_p> & v);

/*------------------------------------------------------------*/
/* assigns a vector to variable "name"                        */
/*------------------------------------------------------------*/
void sage_assign(const Vec<zz_p> & v, const string & name);

/*------------------------------------------------------------*/
/* prints a poly with indeterminate var                       */
/*------------------------------------------------------------*/
void sage_output(const zz_pX & v, const string & var);

/*------------------------------------------------------------*/
/* prints a poly with indeterminate x, as a cast into U       */
/*------------------------------------------------------------*/
void sage_output(const zz_pX & v);

/*------------------------------------------------------------*/
/* assign a poly with indeterminate var to variable "name"    */
/*------------------------------------------------------------*/
void sage_assign(const zz_pX & v, const string & var, const string & name);

/*------------------------------------------------------------*/
/* assign a poly with indeterminate x to variable "name"      */
/*------------------------------------------------------------*/
void sage_assign(const zz_pX & v, const string & name);



#endif
