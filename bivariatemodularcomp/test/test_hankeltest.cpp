#include "mosaic_hankel.h"
#include <sstream>
#include "hermite_pade.h"
#include <NTL/ZZX.h>
using namespace NTL;
using namespace std;

int main(){
  zz_p::init(13);
  Vec<zz_p> v;
  cout << "Enter the function" << endl;
  string s;
  getline(cin,s);
  istringstream iss{s};
  long inp;
  while(iss >> inp)
    v.append(zz_p(inp));
  long prec = v.length();
  Vec<long> type;
  cout << "Enter the type" << endl;
  getline(cin,s);
  iss = istringstream(s);
  while (iss >> inp)
    type.append(inp);
  Vec<hankel> vec_h;
  zz_pX f{1};
  zz_pX f_original;
  conv(f_original, v);
  for (int i = 0; i < type.length(); i++){
    Vec<zz_p> v;
    conv(v,f);
    Vec<zz_p> inp_vec;
    v.SetLength(prec,zz_p{0});
    for (int j = 0; j < v.length(); j++)
      inp_vec.append(v[v.length() - 1 - j]);
    for (int j = 0; j < type[i]; j++)
      inp_vec.append(zz_p{0});
    vec_h.append(hankel(inp_vec,prec,type[i]+1));
    f = f * f_original;
  }
  Mat<zz_p> M;
  for (int i = 0; i < type.length(); i++){
    to_dense(M,vec_h[i]);
    cout << M << endl;
  }
  Vec<Vec<hankel>> init;
  init.append(vec_h);
  mosaic_hankel mh (init);
  to_dense(M,mh);
  cout << M << endl;
  ZZX fp;
  conv(fp, f_original);
  cout << "From Hermite-Pade" << endl;
  hermite_pade(fp,type);
}
