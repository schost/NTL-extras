#include<iostream>
#include<fstream>
#include<sstream>
using namespace std;

/*************************************************
* Generates test cases where there are n entries *
* of n with rank defect 1                        *
*************************************************/
int main(){
  for (long nbits = 1; nbits < 1000; nbits++){
    for (long n = 1; n < 300; n++){
      ostringstream oss;
      oss << "test" << nbits << "_" << n << ".data";
      string filename = oss.str();
      cout << "creating: " << filename << endl;
      ofstream ofs{filename};
      ofs << nbits << endl;
      ofs << 1 << endl;
      for (long i = 0; i < n; i++){
	ofs << n-1 << " ";
      }
    }
  }
}
