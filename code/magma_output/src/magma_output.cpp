#include <NTL/ZZX.h>
#include <NTL/lzz_pX.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>

#include "magma_output.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/* prints a matrix                                            */
/*------------------------------------------------------------*/
void magma_output(const Mat<zz_p> & v){
  if (v.NumRows() == 0){
    cout << "[[]]";
    return;
  }
  cout << "Matrix(GF(" << zz_p::modulus() << "), [";
  for (long i = 0; i < v.NumRows()-1; i++){
    magma_output(v[i]);
    cout << ", ";
  }
  magma_output(v[v.NumRows()-1]);
  cout << "])";
}

/*------------------------------------------------------------*/
/* assigns a vector to variable "name"                        */
/*------------------------------------------------------------*/
void magma_assign(const Mat<zz_p> & v, const string & name){
  cout << name << " := ";
  magma_output(v);
  cout << ";" << endl;
}


/*------------------------------------------------------------*/
/* prints a vector                                            */
/*------------------------------------------------------------*/
void magma_output(const Vec<ZZ> & v){
  if (v.length() == 0){
    cout << "[]";
    return;
  }
  cout << "[";
  for (long i = 0; i < v.length()-1; i++)
    cout << v[i] << ", ";
  cout << v[v.length()-1] << "]";
}

/*------------------------------------------------------------*/
/* assigns a vector to variable "name"                        */
/*------------------------------------------------------------*/
void magma_assign(const Vec<ZZ> & v, const string & name){
  cout << name << " := ";
  magma_output(v);
  cout << ";" << endl;
}

/*------------------------------------------------------------*/
/* initializes U<XX>=Z[XX]                                    */
/*------------------------------------------------------------*/
void magma_init_ZZX(){
  cout << "U<XX> := PolynomialRing(Integers());\n";
}

/*------------------------------------------------------------*/
/* initializes U<XX>=Z[XX]                                    */
/*------------------------------------------------------------*/
void magma_init_QQX(){
  cout << "U<XX> := PolynomialRing(Rationals());\n";
}

/*------------------------------------------------------------*/
/* prints a poly with indeterminate var                       */
/*------------------------------------------------------------*/
void magma_output(const ZZX & v, const string & var){
  cout << "(Parent(" << var << ")!(0)";
  for (long i = 0; i <= deg(v); i++)
    cout << "+(" << coeff(v, i) << ")*" << var << "^" << i;
  cout << ")";
}

/*------------------------------------------------------------*/
/* prints a poly with indeterminate XX                        */
/*------------------------------------------------------------*/
void magma_output(const ZZX & v){
  magma_output(v, "XX");
}
  
/*------------------------------------------------------------*/
/* assign a poly with indeterminate var to variable "name"    */
/*------------------------------------------------------------*/
void magma_assign(const ZZX & v, const string & var, const string & name){
  cout << name << " := ";
  magma_output(v, var);
  cout << ";" << endl;
} 

/*------------------------------------------------------------*/
/* assign a poly with indeterminate XX to variable "name"     */
/*------------------------------------------------------------*/
void magma_assign(const ZZX & v, const string & name){
  magma_assign(v, "XX", name);
}

/*------------------------------------------------------------*/
/* prints a vector                                            */
/*------------------------------------------------------------*/
void magma_output(const Vec<long> & v){
  if (v.length() == 0){
    cout << "[]";
    return;
  }
  cout << "[";
  for (long i = 0; i < v.length()-1; i++)
    cout << v[i] << ", ";
  cout << v[v.length()-1] << "]";
}

/*------------------------------------------------------------*/
/* assigns a vector to variable "name"                        */
/*------------------------------------------------------------*/
void magma_assign(const Vec<long> & v, const string & name){
  cout << name << " := ";
  magma_output(v);
  cout << ";" << endl;
}

/*------------------------------------------------------------*/
/* prints a vector                                            */
/*------------------------------------------------------------*/
void magma_output(const Vec<unsigned long> & v){
  if (v.length() == 0){
    cout << "[]";
    return;
  }
  cout << "[";
  for (long i = 0; i < v.length()-1; i++)
    cout << v[i] << ", ";
  cout << v[v.length()-1] << "]";
}

/*------------------------------------------------------------*/
/* assigns a vector to variable "name"                        */
/*------------------------------------------------------------*/
void magma_assign(const Vec<unsigned long> & v, const string & name){
  cout << name << " := ";
  magma_output(v);
  cout << ";" << endl;
}

/*------------------------------------------------------------*/
/* initializes GF(p)                                          */
/*------------------------------------------------------------*/
void magma_init(){
  cout << "k := GF(" << zz_p::modulus() << ");\n";
}

/*------------------------------------------------------------*/
/* initializes U<x>=GF(p)[x]                                  */
/*------------------------------------------------------------*/
void magma_init_X(){
  cout << "U<x> := PolynomialRing(k);\n";
}

/*------------------------------------------------------------*/
/* prints a vector                                            */
/*------------------------------------------------------------*/
void magma_output(const Vec<zz_p> & v){
  if (v.length() == 0){
    cout << "[]";
    return;
  }
  cout << "[GF(" << zz_p::modulus() << ")|";
  for (long i = 0; i < v.length()-1; i++)
    cout << v[i] << ", ";
  cout << v[v.length()-1] << "]";
}

/*------------------------------------------------------------*/
/* assigns a vector to variable "name"                        */
/*------------------------------------------------------------*/
void magma_assign(const Vec<zz_p> & v, const string & name){
  cout << name << " := ";
  magma_output(v);
  cout << ";" << endl;
}

/*------------------------------------------------------------*/
/* prints a poly with indeterminate var                       */
/*------------------------------------------------------------*/
void magma_output(const zz_pX & v, const string & var){
  cout << "(Parent(" << var << ")!(0)";
  for (long i = 0; i <= deg(v); i++)
    cout << "+(" << coeff(v, i) << ")*" << var << "^" << i;
  cout << ")";
}

/*------------------------------------------------------------*/
/* prints a poly with indeterminate x, as a cast into U       */
/*------------------------------------------------------------*/
void magma_output(const zz_pX & v){
  cout << "U!(";
  magma_output(v.rep);
  cout << ")";
}
  
/*------------------------------------------------------------*/
/* assign a poly with indeterminate var to variable "name"    */
/*------------------------------------------------------------*/
void magma_assign(const zz_pX & v, const string & var, const string & name){
  cout << name << " := ";
  magma_output(v, var);
  cout << ";" << endl;
} 

/*------------------------------------------------------------*/
/* assign a poly with indeterminate x to variable "name"      */
/*------------------------------------------------------------*/
void magma_assign(const zz_pX & v, const string & name){
  magma_assign(v, "x", name);
}

/*------------------------------------------------------------*/
/* prints a vector of polys with indeterminate var            */
/*------------------------------------------------------------*/
void magma_output(const Vec<zz_pX> & v, const string & var){
  cout << "[";
  for (long i = 0; i < v.length(); i++){
    magma_output(v[i]);
    if (i < v.length()-1)
      cout << ", ";
  }
  cout << "]";
}

/*------------------------------------------------------------*/
/* prints a vector of polys with indeterminate x              */
/*------------------------------------------------------------*/
void magma_output(const Vec<zz_pX> & v){
  magma_output(v, "x");
}

/*------------------------------------------------------------*/
/* assign a vector of polys with indet var to variable "name" */
/*------------------------------------------------------------*/
void magma_assign(const Vec<zz_pX> & v, const string & var, const string & name){
  cout << name << ":=";
  magma_output(v, var);
  cout << ";\n";
}

/*------------------------------------------------------------*/
/* assign a vector of polys with indet x to variable "name"   */
/*------------------------------------------------------------*/
void magma_assign(const Vec<zz_pX> & v, const string & name){
  magma_assign(v, "x", name);
}
