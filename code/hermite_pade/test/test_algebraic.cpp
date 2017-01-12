#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <NTL/ZZX.h>

#include "magma_output.h"
#include "lzz_pX_mosaic_hankel.h"
#include "ZZX_extra.h"
#include "ZZXY.h"
#include "ZZ_hermite_pade.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* creates a polynomial F                                     */
/*  deg(F,x)=deg(F,y)=d, ht(F)=b                              */
/* returns the solution with d^2 terms                        */
/*------------------------------------------------------------*/
void random_F_and_solution(ZZXY & F, ZZX & g, ZZ & den_g, long d, long b){
  random(F, b, d, d);

  ZZX F0;
  random(F0, b, F.degY()-1);
  ZZX factor;
  ZZ root{22};
  SetCoeff(factor, 1, 1);
  SetCoeff(factor, 0, -root);
  F0 = F0 * factor;

  for (long i = 0; i <= F.degY(); i++)
    F.coeffX[i] = (F.coeffX[i] << 1) + coeff(F0, i);
    
  F.series_solution(g, den_g, root, 4*d*d);
}

/*------------------------------------------------------------*/
/* reads a polynomial f and an integer vector                 */
/*------------------------------------------------------------*/
void read_f_type(ZZX &F, Vec<long>& type, const string& name){
  string line_f, line_type;
  ifstream file;
  file.open(name);
  long inp;
  Vec<ZZ> f;

  if (! getline(file, line_f)){
    cout << "first line of file should be coefficients of f\n";
    exit(-1);
  }
  istringstream iss_f{line_f};
  while(iss_f >> inp)
    f.append(ZZ(inp));
  F = conv<ZZX>(f);
    
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
/* to test a data file: prog_name filename                    */
/*------------------------------------------------------------*/
int main(int argc, char **argv){

  ZZX f;
  ZZ denom;
  Vec<long> type;
	long mode = 0;

  if (argc == 1){
    ZZXY F;
    long d = 10;
    long b = 2;
    random_F_and_solution(F, f, denom, d, b);
    // magma_init_bi_QQ();
    // magma_assign(F, "F");
    type.SetLength(d+2);
    for (long i = 0; i < d+2; i++)
      type[i] = d+2;
  }
  else {
    read_f_type(f, type, argv[1]);
    if (argc == 3)
    	mode = stoi(argv[2]);
    denom = to_ZZ("4398573498573985749857348957349875349857348957349857349");
  }    
  
  hermite_pade_algebraic hp(f, denom, type, deg(f)+1);
  hp.switch_mode(mode);
  cout << "old type:= " << type << "\n";
  cout << "old dimensions:= " << hp.NumRows() << " " << hp.NumCols() << endl;
  cout << "old rank:= " << hp.Rank() << "\n";
  cout << "old nullity:= " << hp.NumCols() - hp.Rank() << "\n";
  
  if (hp.NumCols() - hp.Rank() > 1){
    Vec<long> new_type = hp.find_new_type();
    hp.init(new_type);
    cout << "new type " << new_type << endl;
    cout << "new dimensions:= " << hp.NumRows() << " " << hp.NumCols() << endl;
    cout << "new rank " << hp.Rank() << endl;
    cout << "new nullity " << hp.NumCols() - hp.Rank() << endl;
  }
  
  Vec<zz_pX> sol_zz_p;
  hp.random_solution_mod_p(sol_zz_p);
  cout << "sol zz_p: " << sol_zz_p << endl;
  
  Vec<ZZX> sol;
  hp.random_solution(sol);
  cout << "sol: " << sol << endl;
}
