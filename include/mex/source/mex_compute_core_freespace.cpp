#include <mex.h>
#include <math.h>

int sign(double x)
{
	if(x >= 0) 
		return 1;
	else
		return -1;
}

int min(int a, int b)
{
	if(a >= b) return b;
	else return a;
}

int max(int a, int b)
{
	if(a >= b) return a;
	else return b;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  
    double *fs = mxGetPr(prhs[0]); 
	double *ifsmask = mxGetPr(prhs[1]);
	double *iptmask = mxGetPr(prhs[2]);
	int minx = (int)*mxGetPr(prhs[3]); 
	int maxx = (int)*mxGetPr(prhs[4]); 
	int miny = (int)*mxGetPr(prhs[5]); 
	int maxy = (int)*mxGetPr(prhs[6]); 

	int h = (int)mxGetDimensions(prhs[0])[0];
	int w = (int)mxGetDimensions(prhs[0])[1]; 

    nlhs = 1;
    plhs[0] = mxCreateDoubleMatrix(4, 1, mxREAL);
	double *X = mxGetPr(plhs[0]);

	double maxcost = 0;
	double maxi, maxj, maxii, maxjj;

	// computing the core-freespace-evidence
	for(int i=minx;i<maxx-1;i++)
	{
		for(int j=miny;j<maxy-1;j++)
		{
			if(fs[i*h+j]>0)
			{
				for(int ii=i+1;ii<maxx;++ii)
				{
					for(int jj=j+1;jj<maxy;++jj)
					{
						
						if(fs[ii*h+jj]>0 && fs[i*h+jj]>0 && fs[ii*h+j]>0)
						{
							
							int ii_ = ii+1;
							int jj_ = jj+1;
							int sum_fs = ifsmask[i*(h+1)+j] + ifsmask[ii_*(h+1)+jj_] - ifsmask[ii_*(h+1)+j] - ifsmask[i*(h+1)+jj_];
							
							if(sum_fs == 0)
							{
								double cost = min(ii-i, jj-j);
								double cost2 = max(ii-i, jj-j);
								cost = cost + 0.001*cost2;;
								if(cost >= maxcost)
								{
									maxcost = cost;
									maxi = i;
									maxj = j;
									maxii = ii;
									maxjj = jj;									
								}
							}
						}
					}

					//}
				}
			}
		}
	}
    
	*(X+0) = maxi+1;
	*(X+1) = maxj+1;
	*(X+2) = maxii+1;
	*(X+3) = maxjj+1;
	
}
