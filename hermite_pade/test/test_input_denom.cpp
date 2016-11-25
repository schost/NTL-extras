#include <iostream>
#include <sstream>
#include "hermite_pade_exact.h"
using namespace std;

int main(){
  ZZX f;
  Vec<long> type;
  string s;
  getline(cin,s);
  istringstream iss(s);
  long x;
  long i = 0;
  cout << "enter f (in one line)" << endl;
  while(iss>>x){
    SetCoeff(f,i++,x);
  }
  cout << "f: " << f << endl;
  getline(cin,s);
  iss = istringstream(s);
  cout << "enter type" << endl;
  while (iss >> x){
    type.append(x);
  }
  cout << "type: " << type << endl;
  ZZ z;
  cout << "enter the denominator" << endl;
  cin >> z;
  cout << "denom: " << z << endl;
  hermite_pade_exact hp(f,z,type,deg(f)+1,10);
  Vec<Vec<ZZ>> sol;
  cout << "new type: " << hp.find_sol(sol) << endl;
  cout << "sol: " << sol << endl;
}
