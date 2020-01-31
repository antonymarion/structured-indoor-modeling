#include <mex.h>
#include <stdlib.h>
#include <cmath>
#include <stdio.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

vector<double> crossproduct(const vector<double> n1, const vector<double> n2) {
      
	vector<double> vec(3);
	vec[0] = n1[1] * n2[2] - n1[2] * n2[1];
    vec[1] = n1[2] * n2[0] - n1[0] * n2[2];
    vec[2] = n1[0] * n2[1] - n1[1] * n2[0];
      return vec;
}

int sign(double x)
{
	if(x >= 0) 
		return 1;
	else
		return -1;
}

const double PI = 3.14159265;
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  
    double *N = mxGetPr(prhs[0]); 
	double *n1 = mxGetPr(prhs[1]); 
	
	int m = (int)mxGetDimensions(prhs[0])[1];
    
    nlhs = 1;
    plhs[0] = mxCreateDoubleMatrix(6, 1, mxREAL);
	double *M = mxGetPr(plhs[0]);
	
	int id;
	double x,y,z;
	double param[4];
	double maxInliers = 0;
	double K = 1.0e5;
	int k=0;
	vector<double> inlierID;
    vector<double> maxinlierID;
	int randomIndex;
	double theta;
	double phi;
	vector<double> nz(3);
	nz[0] = n1[0];
	nz[1] = n1[1];
	nz[2] = n1[2];
	vector<double> nest(3);
	vector<double> nest2(3);
	double ndata[3];
	//mexPrintf("%f %f %f\n",N[0],N[1],N[2]);

	while(k < K){
		inlierID.clear();

		k++;
		// random sampling
		randomIndex = rand()%m;
		// compute the model parameters
		nest[0] = N[3*randomIndex+0];
		nest[1] = N[3*randomIndex+1];
		nest[2] = N[3*randomIndex+2];	

		nest2 =crossproduct(nest,nz);

		// counting inliers and outliers
		double numInliers = 0;
		double ave = 0;
		double std = 0;
		vector<double> err(m);
		for(int i=0;i<m;i++){
			ndata[0] = N[3*i+0];
			ndata[1] = N[3*i+1];
			ndata[2] = N[3*i+2];

			err[i] = min(acos(abs(nest[0]*ndata[0] + nest[1]*ndata[1] + nest[2]*ndata[2])),acos(abs(nest2[0]*ndata[0] + nest2[1]*ndata[1] + nest2[2]*ndata[2])));
			ave = ave + err[i];
		}
		ave = ave/m;

		for(int i=0;i<m;i++){
			std += (err[i]-ave)*(err[i]-ave);
		}
		std = std/m;

		for(int i=0;i<m;i++){
			if(err[i] < 0.01){
				numInliers++;
				inlierID.push_back(i);
			}
		}

		double w = (numInliers-3)/m;


		if(numInliers > maxInliers){
			maxInliers = numInliers;
			*(M+0) = nest[0];
			*(M+1) = nest[1];
			*(M+2) = nest[2];
			*(M+3) = nest2[0];
			*(M+4) = nest2[1];
			*(M+5) = nest2[2];

			double p = max(0.001,pow(w,3));
			K = log(1-0.999)/log(1-p);
			maxinlierID.resize(inlierID.size());
			copy(inlierID.begin(), inlierID.end(), maxinlierID.begin());
		}
		if(k>100000) break;
	}
	
}
