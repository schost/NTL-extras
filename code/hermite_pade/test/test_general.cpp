#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <NTL/ZZX.h>

#include "lzz_pX_mosaic_hankel.h"
#include "ZZ_hermite_pade.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* reads a vecor of polynomials F and an integer vector type  */
/*------------------------------------------------------------*/
void read_vecf_type(Vec<ZZX> &vec_F, Vec<long>& type, const string& name){
  string line_f, line_type;
  ifstream file;
  file.open(name);
  long inp;

  long nb;
  if (! getline(file, line_f)){
    cout << "first line of file should be number of polynomials\n";
    exit(-1);
  }
  istringstream iss_nb{line_f};
  iss_nb >> nb;
  vec_F.SetLength(nb);

  for (long i = 0; i < nb; i++){
    Vec<ZZ> f;
    if (! getline(file, line_f)){
      cout << "not enough polynomials\n";
      exit(-1);
    }
    istringstream iss_f{line_f};
    while(iss_f >> inp)
      f.append(ZZ(inp));
    vec_F[i] = conv<ZZX>(f);
  }

  if (! getline(file, line_type)){
    cout << "last line of file should be type\n";
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

  Vec<ZZX> f;
  Vec<long> type;
  read_vecf_type(f, type, argv[1]);
  cout << f << endl;
  cout << type << endl;

  long s = 0;
  for (long i = 0; i < f.length(); i++)
    s = max(s, deg(f[i]));

  hermite_pade_general hp(f, type, s+1);

  Vec<zz_pX> sol_zz_p;
  hp.random_solution_mod_p(sol_zz_p);
  cout << "sol zz_p: " << sol_zz_p << endl;
  zz_pX res;
  for (long i = 0; i < type.length(); i++)
    res += conv<zz_pX>(f[i]) * sol_zz_p[i];
  cout << "check_zz_p: " << trunc(res, s+1) << endl;
  cout << "k:=GF(" << zz_p::modulus() << ");\n";

  Vec<ZZX> sol;
  hp.random_solution(sol);
  cout << "sol ZZ: " << sol << endl;
  ZZX resZZ;
  for (long i = 0; i < type.length(); i++)
    resZZ += f[i] * sol[i];
  cout << "check ZZ: " << trunc(resZZ, s+1) << endl;

}
