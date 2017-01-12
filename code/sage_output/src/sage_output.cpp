#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

NTL_CLIENT

/*------------------------------------------------------------*/
/* initializes GF(p)                                          */
/*------------------------------------------------------------*/
void sage_init(){
  cout << "k = GF(" << zz_p::modulus() << ")\n";
}

/*------------------------------------------------------------*/
/* initializes U.<x>=GF(p)[x]                                 */
/*------------------------------------------------------------*/
void sage_init_X(){
  cout << "U.<x> = PolynomialRing(k)\n";
}

/*------------------------------------------------------------*/
/* prints a vector                                            */
/*------------------------------------------------------------*/
void sage_output(const Vec<zz_p> & v){
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
void sage_assign(const Vec<zz_p> & v, const string & name){
  cout << name << " = ";
  sage_output(v);
  cout << endl;
}

/*------------------------------------------------------------*/
/* prints a poly with indeterminate x, as a cast into U       */
/*------------------------------------------------------------*/
void sage_output(const zz_pX & v){
  cout << "U(";
  sage_output(v.rep);
  cout << ")";
}

/*------------------------------------------------------------*/
/* prints a poly with indeterminate var                       */
/*------------------------------------------------------------*/
void sage_output(const zz_pX & v, const string & var){
  cout << "(" << var << ".parent()(0)";
  for (long i = 0; i <= deg(v); i++)
    cout << "+(" << coeff(v, i) << ")*" << var << "^" << i;
  cout << ")";
}
 
/*------------------------------------------------------------*/
/* assign a poly with indeterminate var to variable "name"    */
/*------------------------------------------------------------*/
void sage_assign(const zz_pX & v, const string & var, const string & name){
  cout << name << " = ";
  sage_output(v, var);
  cout << endl;
} 

/*------------------------------------------------------------*/
/* assign a poly with indeterminate x to variable "name"      */
/*------------------------------------------------------------*/
void sage_assign(const zz_pX & v, const string & name){
  sage_assign(v, "x", name);
}
