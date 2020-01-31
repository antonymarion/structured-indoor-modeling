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
	int ymin = (int)*mxGetPr(prhs[3]);
	int xmin = (int)*mxGetPr(prhs[4]);
	int xmax = (int)*mxGetPr(prhs[5]);
	int gNumY = (int)mxGetDimensions(prhs[0])[0];
	int gNumX = (int)mxGetDimensions(prhs[0])[1]; 
	vector<double> circump(gNumY);
	double maxheight = -1;
	double maxcost = -1;
	double minheight = -1;

	
	//mexPrintf("H%d W%d\n", gNumY, gNumX);

	// find ceil_line
	for(int y=ymin;y<gNumY;y++)
	{
		
		for(int x=xmin;x<xmax;x++)
		{
			if(wp[x*gNumY+y] > 0 && fsp[x*gNumY+y] == 0) 
			{
				//mexPrintf("wp %f\n", wp[x*gNumY+y] );
				circump[y] += 1;
			}
		}
	}

	double accum = 0;
	for(int y=gNumY-1;y>=ymin;--y)
	{
		double cost = circump[y]-accum;
		accum += circump[y];
		//for(int t=-margin;t<=margin;t++)
		//{
		//	cost += circump[y+t];
		//}

		if(cost >= maxcost)
		{
			//mexPrintf("maxcost %d\n", maxcost);
			maxcost = cost;
			maxheight = y;
		}
	}

	maxcost = -1;
	// find floor_line
	for(int y=ymin;y>=0;--y)
	{
		
		for(int x=0;x<gNumX;x++)
		{
			if(wp[x*gNumY+y] > 0  && fsp[x*gNumY+y] == 0) 
			{
				//mexPrintf("wp %f\n", wp[x*gNumY+y] );
				circump[y] += 1;
			}
		}
	}

	accum = 0;
	for(int y=0;y<ymin;++y)
	{
		double cost = circump[y]-accum;
		accum += circump[y];

		if(cost >= maxcost)
		{
			//mexPrintf("maxcost %d\n", maxcost);
			maxcost = cost;
			minheight = y;
		}
	}


    nlhs = 2;
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	double *maxh = mxGetPr(plhs[0]);
	*maxh = maxheight;	

    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
	double *minh = mxGetPr(plhs[1]);
	*minh = minheight;	
}
