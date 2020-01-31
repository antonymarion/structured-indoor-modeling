#include <mex.h>
#include <cstdio>
#include <math.h>
#include <iostream>
#include <vector>
#include <stdint.h>
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	double *f1 = mxGetPr(prhs[0]);   
	double *f2 = mxGetPr(prhs[1]); 
	double *sum1 = mxGetPr(prhs[2]);
	double *sum2 = mxGetPr(prhs[3]);


	int N1 = mxGetDimensions(prhs[0])[1];
	int D1 = mxGetDimensions(prhs[0])[0];

	int N2 = mxGetDimensions(prhs[1])[1];
	int D2 = mxGetDimensions(prhs[1])[0];

	if (D1!=D2)
	mexErrMsgIdAndTxt("MyProg:InputString","Size of array is invalid.");
	int D = D1;

	nlhs = 1;
    plhs[0] = mxCreateDoubleMatrix(N1*N2, 1, mxREAL);
	double *dist = mxGetPr(plhs[0]);

	
	/* case 1*/
	
	for(int j=0; j<N2; j++)
	{
		for(int i=0; i<N1; i++){
			double c = 0;
			for (int k = 0; k < D; k++)
			{
				if(f1[i*D+k] > 0 &&  f2[j*D+k] == 0){
				c += f1[i*D+k]/sum1[i];
				}

				if(f1[i*D+k] == 0 &&  f2[j*D+k] > 0){
				c += f2[j*D+k]/sum2[j];
				}
			}
			*(dist + j*N1+i) = c;
		}
	}
	


}

