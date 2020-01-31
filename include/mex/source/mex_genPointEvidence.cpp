#include <mex.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <iostream>
#include <vector>
using namespace std;
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *P = mxGetPr(prhs[0]); 
	double *gNum = mxGetPr(prhs[1]);
	double *bSize = mxGetPr(prhs[2]);
	double *b0 = mxGetPr(prhs[3]);
	int m = (int)mxGetDimensions(prhs[0])[1];
	double x[3];
	unsigned int idx, idy, idz;

    nlhs = 1;
	plhs[0] = mxCreateDoubleMatrix(gNum[0]*gNum[1]*gNum[2], 1, mxREAL);
    double *point_evidence =  mxGetPr(plhs[0]);

	for(int i=0; i<m; i++)
	{
		x[0] = P[3*i+0];
		x[1] = P[3*i+1];
		x[2] = P[3*i+2];
		idx = floor(gNum[0]*(x[0]-b0[0])/bSize[0]);
		idy = floor(gNum[1]*(x[1]-b0[1])/bSize[1]);
		idz = floor(gNum[2]*(x[2]-b0[2])/bSize[2]);
		unsigned long long index = idz*gNum[0]*gNum[1]+idx*gNum[1]+idy;
		*(point_evidence+index) = *(point_evidence+index)+1;
	}  	
}
