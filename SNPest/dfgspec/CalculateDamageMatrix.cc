#include <math.h>
#include <sstream>
#include <cctype>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <algorithm>
#include <sstream>
#include <vector>
#include <ctime>

using std::flush;
using std::cout;
using std::endl;
using std::cerr;
using std::cin;
using std::ifstream;
using std::istream;
using std::ofstream;
using std::stringstream;
using std::ios;
using std::string;
using std::vector;

// Main function. Read in sequences. Remember names.
int main(int argc, char *argv[]){
  if(argc!=3){
    cerr<<"Usage: "<<argv[0]<<" tau delta\n";
    return 1;
  }
  double tau=atof(argv[1]);
  double delta=atof(argv[2]);
  
  if(tau<0 || tau>=1 || delta<0 || delta>=1){
	cerr<<"Illegal parameter settings. Use 0<=tau<1 and 0<=delta<1:\ntau="<<tau<<"\ndelta="<<delta<<endl;
      return 1;
    }
												      
  double errormatrix[4][4]={{1.0-tau-delta,3.0/7.0*tau,3.0/7.0*tau,1.0/7.0*tau},
			    {3.0/7.0*tau,1.0-tau,1.0/7.0*tau,3.0/7.0*tau+delta},
			    {3.0/7.0*tau+delta,1.0/7.0*tau,1.0-tau,3.0/7.0*tau},
			    {1.0/7.0*tau,3.0/7.0*tau,3.0/7.0*tau,1.0-tau-delta}};
  
  double PHREDMATRIX[4][200];
  double PT,PE;
  char base[4]={'A','C','G','T'};
  for(int i=0;i<4;i++){
    for(int j=0;j<50;j++){
      //      cerr<<base[i]<<(j+1)<<" ";
      PE=pow(10,-0.1*(j+1.0));
      PT=1.0-PE;
      for(int k=0;k<4;k++){
	if(i==k) PHREDMATRIX[k][j+50*i]=PT;
	else PHREDMATRIX[k][j+50*i]=PE/3;
      }

      // Correct for error rates
      double tempvector[4]={PHREDMATRIX[0][j+50*i],PHREDMATRIX[1][j+50*i],PHREDMATRIX[2][j+50*i],PHREDMATRIX[3][j+50*i]};
      for(int k=0;k<4;k++){
	PHREDMATRIX[k][j+50*i]=errormatrix[k][0]*tempvector[0]+errormatrix[k][1]*tempvector[1]+errormatrix[k][2]*tempvector[2]+errormatrix[k][3]*tempvector[3];
      }
    }
  }

      /*  cout<<"Error matrix:\t(";
  for(int i=0;i<4;i++){
    cout<<"(";
    for(int j=0;j<4;j++){
      cout<<errormatrix[i][j]<<", ";
    }
    cout<<"),\n";
    }*/

  cout<<"POT_MAT:\t[4, 200] (";
  for(int i=0;i<4;i++){
    cout<<"(";
    for(int j=0;j<200;j++){
      cout<<PHREDMATRIX[i][j]<<", ";
    }
    cout<<"),\n";
  }
  cout<<")\n";

	/*  cerr<<endl;
  for(int i=1;i<=50;i++) cerr<<"N"<<i<<"=A"<<i<<" C"<<i<<" G"<<i<<" T"<<i<<"; ";
  cerr<<endl;
	*/
  return 0;

}
