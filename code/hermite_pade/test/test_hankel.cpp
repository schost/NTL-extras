#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <NTL/ZZX.h>

#include "lzz_pX_mosaic_hankel.h"
#include "ZZ_hermite_pade.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* reads a polynomial f and an integer vector                 */
/*------------------------------------------------------------*/
void read_f_type(Vec<zz_p> &f, Vec<long>& type, const string& name){
  string line_f, line_type;
  ifstream file;
  file.open(name);
  long inp;

  if (! getline(file, line_f)){
    cout << "first line of file should be coefficients of f\n";
    exit(-1);
  }
  istringstream iss_f{line_f};
  while(iss_f >> inp)
    f.append(zz_p(inp));
    
  if (! getline(file, line_type)){
    cout << "first line of file should be type\n";
    exit(-1);
  }
  istringstream iss_type = istringstream(line_type);
  while (iss_type >> inp)
    type.append(inp);

  file.close();
}

/*------------------------------------------------------------*/
/* usage:                                                     */
/* prog_name filename                                         */
/*------------------------------------------------------------*/
int main(int argc, char **argv){

  if (argc == 1){
    cout << "usage: prog_name filename\n";
    return 0;
  }

  zz_p::init(13);
  Vec<zz_p> v;
  Vec<long> type;
  read_f_type(v, type, argv[1]);
  long prec = v.length();

  Vec<hankel> vec_h;
  zz_pX f{1};
  zz_pX f_original;
  conv(f_original, v);
  for (int i = 0; i < type.length(); i++){
    Vec<zz_p> v;
    conv(v, f);
    Vec<zz_p> inp_vec;
    v.SetLength(prec, zz_p{0});
    for (int j = 0; j < v.length(); j++)
      inp_vec.append(v[v.length() - 1 - j]);
    for (int j = 0; j < type[i]; j++)
      inp_vec.append(zz_p{0});
    vec_h.append(hankel(inp_vec, prec,type[i]+1));
    f = f * f_original;
  }

  Mat<zz_p> M;
  for (int i = 0; i < type.length(); i++){
    to_dense(M, vec_h[i]);
    cout << M << endl;
  }
  Vec<Vec<hankel>> init;
  init.append(vec_h);
  mosaic_hankel mh (init);
  to_dense(M, mh);
  cout << M << endl;
}
