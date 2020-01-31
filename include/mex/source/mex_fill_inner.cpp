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
  
    double *xlist =  mxGetPr(prhs[0]);
	double *ylist =  mxGetPr(prhs[1]);
	int gNumX = (int) *mxGetPr(prhs[2]);
	int gNumY = (int) *mxGetPr(prhs[3]); 
    int pNum = (int) mxGetDimensions(prhs[0])[0];
    nlhs = 1;
    plhs[0] = mxCreateDoubleMatrix(gNumX*gNumY, 1, mxREAL);
	vector<double> mask(gNumX*gNumY);
	vector<double> mask2(gNumX*gNumY);
	double *maskOut = mxGetPr(plhs[0]);

	//mexPrintf("...making mask\n");
	int x0, y0, x1, y1, x, y, x_, y_;
	// generate mask
	for(int i=0;i<pNum-1;i++){
		x0 = xlist[i]-1;
		y0 = ylist[i]-1;
		x1 = xlist[i+1]-1;
		y1 = ylist[i+1]-1;
		double vx = x1-x0;
		double vy = y1-y0;
		double nv = sqrt(vx*vx+vy*vy);

		if(nv>0)
		{
			vx = vx/nv;
			vy = vy/nv;
			double count = 0;
			while(1){
				x = floor(x0 + 0.5*count*vx);
				y = floor(y0 + 0.5*count*vy);
				if(x< 0 || x>= gNumX || y < 0 || y >= gNumY) break;
				//mexPrintf("[%d %d] [%d %d] %d %d (%d)\n", x0, y0, x1, y1, x,y, gNumY);
				mask[x*gNumY+y] = 1;
				count++;
				if(x==x1 && y==y1) break;
			}
		}
	}

	//mexPrintf("...filling the inner object\n");
	////// computing the inner object
	for(int i=0;i<pNum-1;i++){
	x0 = xlist[i]-1;
	y0 = ylist[i]-1;
	x1 = xlist[i+1]-1;
	y1 = ylist[i+1]-1;
	double vx = x1-x0;
	double vy = y1-y0;
	double nv = sqrt(vx*vx+vy*vy);
	vx = vx/nv;
	vy = vy/nv;
	double vx2 = vy;
	double vy2 = -vx;
	int count = 0;
	while(1){

		x = floor(x0 + 0.5*count*vx);
		y = floor(y0 + 0.5*count*vy);
		if(x< 0 || x >= gNumX || y < 0 || y >= gNumY) break;
		mask2[x*gNumY+y] = 1;
		if((x!=x0 || y!=y0) & (x!=x1 || y!=y1)){
		int count2 = 1;
		while(1){				
			x_ = x + 0.5*count2*vx2;
			y_ = y + 0.5*count2*vy2;
			if(x_< 0 || x_ >= gNumX || y_ < 0 || y_ >= gNumY) break;
			if(x_!=x0 && y_!=y0 && mask[x_*gNumY+y_]==1) break;				
			mask2[x_*gNumY+y_] = 1;
			count2++;
		}
		
		}
		

		if(x==x1 && y==y1) break;
		count++;
		
	}
	}
	for(int i=0;i<gNumX*gNumY;i++)
	{
		if(mask[i]==1||mask2[i]==1)
		{
			*(maskOut+i) = 1;
		}
		else
		{
			*(maskOut+i) = 0;
		}
	}
	
}
