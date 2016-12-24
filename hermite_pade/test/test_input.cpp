#include <iostream>
#include <sstream>
#include "hermite_pade_exact.h"
using namespace std;

int main(){
  ZZX f;
  Vec<long> type;
  string s;
  getline(cin,s);
  cout << "line: " << s << endl;
  istringstream iss(s);
  ZZ x;
  long t;
  long i = 0;
  while(true){
  	iss>>x;
  	if (iss.eof()) break;
  	if (iss.fail() && !iss.eof()){
  		iss.clear();
  		iss.ignore();
  	}
  	
    cout << x << endl;
    SetCoeff(f,i++,x);
  }
  cout << "f: " << f << endl;
  getline(cin,s);
  iss = istringstream(s);
  while (iss >> t){
    type.append(t);
  }
  cout << "type: " << type << endl;
  hermite_pade_exact hp(f,type,deg(f)+1,10);
  Vec<Vec<ZZ>> sol;
  cout << "new type: " << hp.find_sol(sol) << endl;
  cout << "sol: " << sol << endl;
}
