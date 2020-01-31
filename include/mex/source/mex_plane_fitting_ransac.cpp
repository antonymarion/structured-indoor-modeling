#include <mex.h>
#include <stdlib.h>
#include <cmath>
#include <stdio.h>
#include <time.h>
#include "Eigen/Dense"
#include <iostream>
#include <vector>
using namespace Eigen;
using namespace std;
//#define pi = 3.141592

double min(double x, double y)
{
	if(x > y)
		return y;
	else
		return x;
	
}

double max(double x, double y)
{
	if(x > y)
		return x;
	else
		return y;
	
}

double norm(double x, double y, double z){

	return sqrt(x*x+y*y+z*z);

}

//double abs(double x)
//{	
//	if(x >= 0) return x;
//	else return -x;	
//}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    const double *X = mxGetPr(prhs[0]);
	const double *Y = mxGetPr(prhs[1]);
    const double *Z = mxGetPr(prhs[2]);
	
    const double *n0 = mxGetPr(prhs[3]);
    const double delta = *mxGetPr(prhs[4]);
	const double t = *mxGetPr(prhs[5]);
	const int m = mxGetDimensions(prhs[0])[0];
	
	nlhs = 2;
	plhs[0] = mxCreateDoubleMatrix(4, 1, mxREAL);
	double *maxParam = mxGetPr(plhs[0]);
	plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
	double *iter = mxGetPr(plhs[2]);

	srand((unsigned int)time(NULL));
	
	int id;
	double x,y,z;
	double param[4];
	double maxInliers = 0;
	double K = 1.0e5;
	int k=0;
	vector<double> inlierID;
    vector<double> maxinlierID;
	while(k < K){
		inlierID.clear();

		k++;
		// random sampling
		int count = 0;
		int idlist[3];
		bool flag;
		while(count<3){
			id = rand()%m;

			// all point should be extracted from different position
			flag = true;
			for(int j=0;j<count;j++){
				if(id==idlist[j]) flag = false;
			}
			if(flag == false) continue;

			x = X[id];
			y = Y[id];
			z = Z[id];
			idlist[count] = id;
			count++;
		}

		// compute the model parameters
		MatrixXd M(3,4);
		M << X[idlist[0]], Y[idlist[0]], Z[idlist[0]], 1, X[idlist[1]], Y[idlist[1]], Z[idlist[1]], 1, X[idlist[2]], Y[idlist[2]], Z[idlist[2]], 1;
		MatrixXd MtM = M.transpose()*M; 
		
		SelfAdjointEigenSolver<MatrixXd> es(MtM);
		MatrixXd V = es.eigenvectors();
		for(int i=0;i<4;i++){
			param[i] = V(i,0);
		}

		// constraint check
		double nx, ny, nz, nm, theta;
		nx = param[0];
		ny = param[1];
		nz = param[2];
		nm = norm(nx,ny,nz);
		nx = nx/nm;
		ny = ny/nm;
		nz = nz/nm;
		theta = acos(nx*n0[0] + ny*n0[1] + nz*n0[2]);
		if(theta > delta) continue;

		// counting inliers and outliers
		double numInliers = 0;
		double ave = 0;
		double std = 0;
		vector<double> err(m);
		for(int i=0;i<m;i++){
			err[i] = (param[0]*X[i] + param[1]*Y[i] + param[2]*Z[i] + param[3])/param[2];
			err[i] = abs(err[i]);
			ave = ave + err[i];
		}
		ave = ave/m;

		for(int i=0;i<m;i++){
			std += (err[i]-ave)*(err[i]-ave);
		}
		std = std/m;

		for(int i=0;i<m;i++){
			if(err[i]*err[i] < t*std){
				numInliers++;
				inlierID.push_back(i);
			}
		}

		double w = (numInliers-3)/m;


		if(numInliers > maxInliers){
			maxInliers = numInliers;
			for(int j=0;j<4;j++) *(maxParam+j) = param[j];
			double p = max(0.001,pow(w,3));
			K = log(1-0.999)/log(1-p);
			maxinlierID.resize(inlierID.size());
			copy(inlierID.begin(), inlierID.end(), maxinlierID.begin());
		}
		if(k>100000) break;
	}

	*iter = k;

	plhs[1] = mxCreateDoubleMatrix(maxinlierID.size(), 1, mxREAL);
	double *IID = mxGetPr(plhs[1]);

	for(int i=0;i<maxinlierID.size();i++){
		*(IID+i) = maxinlierID[i]+1;
	}
	//mexPrintf("%f\n",maxInliers);
}
	
