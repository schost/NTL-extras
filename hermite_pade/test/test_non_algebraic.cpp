#include <iostream>
#include <sstream>
#include "hermite_pade_non_algebraic.h"
using namespace std;
using namespace NTL;

int main(){
  int num_fs;
  Vec<ZZX> fs;
  string s;
  cout << "Enter the number of functions" << endl;
  getline(cin,s);
  num_fs = stoi(s);
  for (int i = 0; i < num_fs; i++){
    cout << "Enter the function in one line" << endl;
    getline(cin,s);
    istringstream iss{s};
    int num;
    Vec<ZZ> f;
    while(iss>>num){
      f.append(ZZ(num));
    }
    fs.append(conv<ZZX>(f));
  }
  cout << "Enter the type in one line" << endl;
  getline(cin,s);
  istringstream iss{s};
  Vec<long> type;
  while (iss >> num_fs)
    type.append(num_fs);
  cout << "fs: " << fs << endl;
  cout << "type: " << type << endl;
  hermite_pade_non_algebraic(fs, type, 3);
}
