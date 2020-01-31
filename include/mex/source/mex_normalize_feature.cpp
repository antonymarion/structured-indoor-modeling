#include <mex.h>
#include <vector>
#include <math.h>
using namespace std;

int sign(double x)
{
	if(x >= 0) 
		return 1;
	else
		return -1;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  
    double *feature = mxGetPr(prhs[0]); // grid of the bounding box
	double *countFeatures = mxGetPr(prhs[1]); // FS evidence  (size: gNum(1)*gNum(2))
	int num_seed = mxGetDimensions(prhs[0])[0];
	int num_target = mxGetDimensions(prhs[0])[1];
    
    nlhs = 1;
    plhs[0] = mxCreateDoubleMatrix(num_seed*num_target, 1, mxREAL);
	double *F = mxGetPr(plhs[0]);

	mexPrintf("numSeed %d numTarget %d\n",num_seed, num_target);
	for(int i=0;i<num_seed;i++){
		double acum_val = 0;
		vector<double> val(num_target);
		for(int j=0;j<num_target;j++){
			double fval = feature[num_seed*j+i];
			double cnt = countFeatures[j];
			if(fval>0){
				acum_val += fval/cnt;
				val[j] = fval/cnt;
			}
			else{
				val[j] = 0;
			}
		}
		for(int j=0;j<num_target;j++){
			*(F+num_seed*j+i) = val[j]/acum_val;
		}
	}
	

	
		
	
}
