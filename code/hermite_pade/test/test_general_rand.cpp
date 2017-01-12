#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <NTL/ZZX.h>
#include <NTL/ZZ.h>

#include "lzz_pX_mosaic_hankel.h"
#include "ZZ_hermite_pade.h"

NTL_CLIENT

void rand_gen(Vec<ZZX> &f, Vec<long> &type, long& nbits, char* filename){
	ifstream in{filename};
	in >> nbits;
	
	string s;
	getline(in, s);
	long n;
	long ncols = 0;
	long blocks = 0;
	while (in >> n){
		type.append(n);
		ncols += n+1;
		blocks++;
	}
		
	// randomly generating the entries
	ZZ entry;
	long nrows = ncols - 10;
	for (long j = 0; j < blocks; j++){
		ZZX y;
		for (long i = 0; i < nrows; i++){
			RandomLen(entry, nbits);
			SetCoeff(y,i,entry);
		}
		f.append(y);
	}
	cout << "f: " << f << endl;
	cout << "type: " << type << endl;
}

void run_DAC(const Vec<ZZX> &f, const Vec<long> &type, long prec){
	cout << "Running DAC" << endl;
	hermite_pade_general hp(f, type, prec);
	Vec<ZZX> sol;
  hp.random_solution(sol);
  //cout << "sol ZZ: " << sol << endl;
  ZZX resZZ;
  for (long i = 0; i < type.length(); i++)
    resZZ += f[i] * sol[i];
  cout << "check ZZ: " << trunc(resZZ, prec) << endl;
}

void run_Dixon(const Vec<ZZX> &f, const Vec<long> &type, long prec){
	cout << "Running Dixon" << endl;
	hermite_pade_general hp(f, type, prec);
	hp.switch_mode(1);
	Vec<ZZX> sol;
  hp.random_solution(sol);
  //cout << "sol ZZ: " << sol << endl;
  ZZX resZZ;
  for (long i = 0; i < type.length(); i++)
    resZZ += f[i] * sol[i];
  cout << "check ZZ: " << trunc(resZZ, prec) << endl;
}

void run_CRT(const Vec<ZZX> &f, const Vec<long> &type, long prec, long nbits){
	long times = nbits / 62 + 1;
	cout << "Running CRT " << times << " times" << endl;
	for (long i = 0; i < times; i++){
		hermite_pade_general hp(f, type, prec,i);

  	Vec<zz_pX> sol_zz_p;
  	hp.random_solution_mod_p(sol_zz_p);
 	  cout << "sol zz_p: " << sol_zz_p << endl;
  	zz_pX res;
  	for (long i = 0; i < type.length(); i++)
    	res += conv<zz_pX>(f[i]) * sol_zz_p[i];
  	cout << "check_zz_p: " << trunc(res, prec) << endl;
  	cout << "k:=GF(" << zz_p::modulus() << ");\n";
	}
}

int main(int argc, char **argv){

  if (argc < 3){
    cout << "usage: prog_name filename mode\n";
    return 0;
  }

	long mode = stoi(argv[2]);
  Vec<ZZX> f;
  Vec<long> type;
  long nbits;
  rand_gen(f, type, nbits, argv[1]);
	
  long s = 0;
  for (long i = 0; i < f.length(); i++)
    s = max(s, deg(f[i]));
    
  double t;
  if (mode == 0){
  	t = GetTime();
  	run_CRT(f,type,s+1,nbits);
  	cout << "Took: " << GetTime()-t << endl;
  }else if(mode == 1){
  	t = GetTime();
  	run_DAC(f,type,s+1);
  	cout << "Took: " << GetTime()-t << endl;
  }else{
  	t = GetTime();
  	run_Dixon(f,type,s+1);
  	cout << "Took: " << GetTime()-t << endl;
  }
  

/*  hermite_pade_general hp(f, type, s+1);

  Vec<zz_pX> sol_zz_p;
  hp.random_solution_mod_p(sol_zz_p);
  cout << "sol zz_p: " << sol_zz_p << endl;
  zz_pX res;
  for (long i = 0; i < type.length(); i++)
    res += conv<zz_pX>(f[i]) * sol_zz_p[i];
  cout << "check_zz_p: " << trunc(res, s+1) << endl;
  cout << "k:=GF(" << zz_p::modulus() << ");\n";
*/	
 
}
