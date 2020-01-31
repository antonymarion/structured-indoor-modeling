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
	double *N = mxGetPr(prhs[1]); 
	double *gNum = mxGetPr(prhs[2]);
	double *leaf_size = mxGetPr(prhs[3]);
	double *b0 = mxGetPr(prhs[4]);
	int m = (int)mxGetDimensions(prhs[0])[1];
	double x[3], n[3];
	unsigned int idx, idy, idz;

    nlhs = 4;
	plhs[0] = mxCreateDoubleMatrix(gNum[0]*gNum[1]*gNum[2], 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(gNum[0]*gNum[1]*gNum[2], 1, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(gNum[0]*gNum[1]*gNum[2], 1, mxREAL);
	plhs[3] = mxCreateDoubleMatrix(gNum[0]*gNum[1]*gNum[2], 1, mxREAL);
    double *point_evidence =  mxGetPr(plhs[0]);
	double *nx =  mxGetPr(plhs[1]);
	double *ny =  mxGetPr(plhs[2]);
	double *nz =  mxGetPr(plhs[3]);

	for(int i=0; i<m; i++)
	{
		x[0] = P[3*i+0];
		x[1] = P[3*i+1];
		x[2] = P[3*i+2];
		n[0] = N[3*i+0];
		n[1] = N[3*i+1];
		n[2] = N[3*i+2];
		idx = floor((x[0]-b0[0])/leaf_size[0]);
		idy = floor((x[1]-b0[1])/leaf_size[1]);
		idz = floor((x[2]-b0[2])/leaf_size[2]);
		unsigned long long index = idz*gNum[0]*gNum[1]+idx*gNum[1]+idy;
		*(point_evidence+index) = *(point_evidence+index)+1;
		
		// compute the average vector
		*(nx+index) = n[0] + *(nx+index);
		*(ny+index) = n[1] + *(ny+index);
		*(nz+index) = n[2] + *(nz+index);
	}

	for(int i=0; i<gNum[0]*gNum[1]*gNum[2]; i++)
	{
		if(*(point_evidence+i)>0)
		{
			n[0] = *(nx+i);
			n[1] = *(ny+i);
			n[2] = *(nz+i);
			double nm = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
			*(nx+i) = n[0]/nm;
			*(ny+i) = n[1]/nm;
			*(nz+i) = n[2]/nm;
		}
	}
}
