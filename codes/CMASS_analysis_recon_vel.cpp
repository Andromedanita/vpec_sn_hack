#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include "spline.h"
using namespace std;
#define PI 3.14159265

class star {
public:
  double ra;
  double dec;
  double r;
  double z;
  double v;
  double getTheta(){
    return ((90 - dec)* PI)/180;
  }
  double getDec(){
    return (dec* PI)/180;
  }
  double getPhi(){
    return (ra * PI)/180;
  }
};


double separation(star a, star b);
double cij(star a, star b);
double vi_vj_diff(star a, star b);
double cos_angle(star a, star b);

int main( int argc, char* argv[]){
  // Raising error if number of input files is incorrect
  if (argc !=2 ) {
    cout << "Please provide the file name" << endl;
    exit(1);
  }
  
  string fileName(argv[1]);  // file including (RA, DEC, z, r, vtot) values
  
  cout << "File name is : " << fileName << endl;
  
  ifstream fin(fileName.c_str());
  ofstream fout("output.txt");
  int      numStars = 0;
  
  fin  >> numStars;
  cout << "Number of Stars = " << numStars << endl;
  star *stars = new star[numStars];
  for (int i =0 ; i < numStars; i++){
    float temp;
    float tempz;
    float tempmu;
    float tempnoise;
    // setting values in the files to the star class
    fin  >> stars[i].ra >> stars[i].dec >> temp >> stars[i].r >> stars[i].v;
  }
  
  // computation
  const int num_bins = 20;  // number of separation bins
  int       max_dist = 100; // maximum separation (100 Mpc in this case)
  double    bin_size = (double)(max_dist)/num_bins;
  
  double numVals[num_bins];
  double denomVals[num_bins];
  double num_vals_in_each_bin[num_bins];
  cout << "start computation" << endl;
    
  cij(stars[0], stars[1]);  
    
  int printInterval = 500; // for printing progress every 1000 computations

  for (int i = 0 ; i < numStars ; i++){
    if ( i% printInterval == 0)
      cout << "Iteration " << i << endl;
    for (int j = i+1 ; j < numStars ; j++) {
      double sep = separation(stars[i],stars[j]); // calculating separation

      // Doing calculation of cij and vij difference for separtions smaller than 100 Mpc
      if (sep < 100.0){
	int curbin         = (sep/bin_size); // number of separation bin each pair belongs to

	//cout << "Current bin before modification is" << curbin << endl;
	if (curbin == num_bins){
			 curbin -=1;
			 cout << "changing bin number!" << endl;
	  }

	if (curbin>=num_bins){
	  cout << "error" << endl;
	  cout << sep << " " << bin_size << " " << curbin << endl; 
	}
	
	double num_val     = vi_vj_diff(stars[i], stars[j]) * cij(stars[i], stars[j]); // (vij_dff * cij)
	double denom_val   = pow(cij(stars[i], stars[j]),2);                           // (cij^2)
	numVals[curbin]   += num_val;
	denomVals[curbin] += denom_val;
	num_vals_in_each_bin[curbin] += 1;
	fout << i << "\t\t" << j << "\t\t" << sep << "\t\t" << cij(stars[i], stars[j]) << "\t\t" << cos_angle(stars[i], stars[j]) << "\t\t" << curbin << endl;
      }
    }
  }
  ofstream fout2("vpec_Cmass.txt");
  for (int i = 0; i < num_bins; i++){
    double rval = ((i*bin_size) + ((i+1)*bin_size))/2;
    double fraction  = (numVals[i] / denomVals[i]); 
    fout2 << rval << "\t" << fraction << "\t" << num_vals_in_each_bin[i] << endl;
    cout << "number in each bin is:" << num_vals_in_each_bin[i] << endl;
    }

  fout2.close();
  
  fout.close();
  fin.close();
  return 0;
}

double separation(star a, star b){
  double thetaA = a.getTheta();
  double phiA   = a.getPhi();
  double thetaB = b.getTheta();
  double phiB   = b.getPhi();
  double angles = sin(thetaA) * sin(thetaB) * cos(phiA - phiB) + cos(thetaA) *cos(thetaB);
  double d      = sqrt( a.r * a.r + b.r * b.r - (2 * a.r * b.r * angles));
  return d;
}

double cij(star a, star b){
  double decA  = a.getDec();
  double decB  = b.getDec();
  double raA   = a.getPhi();
  double raB   = b.getPhi();
  double cosD  = sin(decA) * cos(raA) * sin(decB) * cos(raB) + sin(decA) * sin(decB) * sin(raA) * sin(raB)  + cos(decA) * cos(decB); 

  double num = (a.r - b.r) * (1 + cosD);
  double denom = 2 * sqrt( (a.r * a.r) + (b.r * b.r) - (2 * a.r * b.r * cosD));
  
  return num/denom ;
}


double cos_angle(star a, star b){
  double decA  = a.getDec();
  double decB  = b.getDec();
  double raA   = a.getPhi();
  double raB   = b.getPhi();
  double cosD  = sin(decA) * cos(raA) * sin(decB) * cos(raB) + sin(decA) * sin(decB) * sin(raA) * sin(raB)  + cos(decA) * cos(decB);

  return cosD;
}


double vi_vj_diff(star a, star b){
  double vi   = a.v;
  double vj   = b.v;
  double diff = vi - vj;
  return diff;
  }
