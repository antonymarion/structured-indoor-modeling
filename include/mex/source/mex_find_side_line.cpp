#include <mex.h>
#include <math.h>
#include <vector>
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
  
    double *wp = mxGetPr(prhs[0]); 
	double *fsp = mxGetPr(prhs[1]); 
	int margin = (int)*mxGetPr(prhs[2]);
	int xmin = (int)*mxGetPr(prhs[3]);
	int ymin = (int)*mxGetPr(prhs[4]);
	int ymax = (int)*mxGetPr(prhs[5]);
	int gNumY = (int)mxGetDimensions(prhs[0])[0];
	int gNumX = (int)mxGetDimensions(prhs[0])[1]; 
	vector<double> circump(gNumX);
	double maxright = -1;
	double maxcost = -1;
	double minleft = -1;

	
	//mexPrintf("H%d W%d\n", gNumY, gNumX);

	// find right_line
	for(int x=0;x<gNumX;x++)
	{
		
		for(int y=ymin;y<ymax;y++)
		{
			if(wp[x*gNumY+y] > 0 && fsp[x*gNumY+y] == 0) 
			{
				//mexPrintf("wp %f\n", wp[x*gNumY+y] );
				circump[x] += 1;
			}
		}
	}

	double accum = 0;
	for(int x=gNumX-1;x>=0;--x)
	{
		double cost = circump[x]-accum;
		accum += circump[x];
		//for(int t=-margin;t<=margin;t++)
		//{
		//	cost += circump[y+t];
		//}

		if(cost >= maxcost)
		{
			//mexPrintf("maxcost %d\n", maxcost);
			maxcost = cost;
			maxright = x;
		}
	}

	maxcost = -1;
	// find floor_line
	for(int x=gNumX-1;x>=0;--x)
	{
		
		for(int y=ymin;y<ymax;y++)
		{
			if(wp[x*gNumY+y] > 0  && fsp[x*gNumY+y] == 0) 
			{
				//mexPrintf("wp %f\n", wp[x*gNumY+y] );
				circump[x] += 1;
			}
		}
	}

	accum = 0;
	for(int x=0;x<gNumX;++x)
	{
		double cost = circump[x]-accum;
		accum += circump[x];

		if(cost >= maxcost)
		{
			//mexPrintf("maxcost %d\n", maxcost);
			maxcost = cost;
			minleft = x;
		}
	}


    nlhs = 2;
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	double *maxh = mxGetPr(plhs[0]);
	*maxh = maxright;	

    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
	double *minh = mxGetPr(plhs[1]);
	*minh = minleft;	
}
