#include <mex.h>
#include <math.h>
#include <iostream>
#include <queue>
#include <vector>
using namespace std;


/*******************
the direction of path is counter-clockwise path up : core : left
**********************/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	const double *core_mask = mxGetPr(prhs[0]); // core free-space evidence
	const double *points = mxGetPr(prhs[1]); // candidate centers of start-end line
	const double *path_evidence = mxGetPr(prhs[2]); // if the first shortest path result provided, each point is weighted by the cost 
	const double *point_evidence = mxGetPr(prhs[3]); // point evidence
	const double *invalid = mxGetPr(prhs[4]);

	const int h = mxGetDimensions(prhs[0])[0];
    const int w = mxGetDimensions(prhs[0])[1];
	const int numpoints = mxGetDimensions(prhs[1])[1];

	const int dx[] = {-1, 0, 1, 0};
	const int dy[] = {0, 1, 0, -1};
	
	const int sx[] = {0, 1, 0, -1};
	const int sy[] = {1, 0, -1, 0};
	const int ex[] = {0, -1, 0, 1};
	const int ey[] = {-1, 0, 1, 0};

	nlhs = 3;
    plhs[0] = mxCreateDoubleMatrix(2, 1, mxREAL);
	double *xStart = mxGetPr(plhs[0]);
	
	plhs[1] = mxCreateDoubleMatrix(2, 1, mxREAL);
	double *xGoal = mxGetPr(plhs[1]);

	plhs[2] = mxCreateDoubleMatrix(h*w, 1, mxREAL);
	double *mask = mxGetPr(plhs[2]);
	for(int i=0;i<h*w;i++) *(mask+i) = core_mask[i]; 


	
	bool flag1, flag2, flag3;
	double maxcost = -1;
	int max_d;
	int xc[2];
	int xs[2];
	int xe[2];
	int x[2];
	int max_xc[2];
	int max_xs[2];
	int max_xe[2];
	int p;
	double cost;
	double cost2;
	for(int i=0;i<numpoints;i++)
	{
		xc[0] = points[2*i+0];
		xc[1] = points[2*i+1];

		for(int d=0;d<4;d++)
		{
			xs[0] = xc[0] + sx[d];
			xs[1] = xc[1] + sy[d];
			xe[0] = xc[0] + ex[d];
			xe[1] = xc[1] + ey[d];
			flag1 = false;
			flag2 = true;

			p = 0;
			cost2 = 0;
			while(1) // in room
			{
				x[0] = xc[0] + p*dx[d];
				x[1] = xc[1] + p*dy[d];

				
				if(x[0]>=w || x[0]<0 || x[1]>=h || x[1]<0)
				{
					flag1 = false;
					break;						
				}

				if(invalid[x[0]*h+x[1]]>0)
				{
					flag1 = false;
					break;
				}				

				if(*(mask+x[0]*h+x[1])>0)
				{
					flag1 = true;
					break;
				}
				p++;
			}


			p = 0;
			while(1) //out room
			{
				x[0] = xc[0] - p*dx[d];
				x[1] = xc[1] - p*dy[d];



				if(x[0]>=w || x[0]<0 || x[1]>=h || x[1]<0)
				{
					break;						
				}

				if(*(mask+x[0]*h+x[1])>0)
				{
					flag2 = false;
					break;
				}

				//if(invalid[x[0]*h+x[1]]>0)
				//{
				//	flag2 = false;
				//	break;
				//}	
				
				cost2+= point_evidence[x[0]*h+x[1]];

				p++;
			}

			if(flag1==true && flag2==true)
			{
				cost = point_evidence[xc[0]*h+xc[1]] + point_evidence[xs[0]*h+xs[1]] + point_evidence[xe[0]*h+xe[1]] + path_evidence[i] - cost2;

				if(cost > maxcost)
				{
					maxcost=cost;
					max_d = d;
					max_xc[0] = xc[0];
					max_xc[1] = xc[1];
					max_xs[0] = xs[0];
					max_xs[1] = xs[1];
					max_xe[0] = xe[0];
					max_xe[1] = xe[1];
				}
			}

		}
	}
	
	*(xStart+0) = max_xs[0];
	*(xStart+1) = max_xs[1];
	*(xGoal+0) = max_xe[0];
	*(xGoal+1) = max_xe[1];

	p=0;
	while(1){
		x[0] = max_xc[0] + p*dx[max_d];
		x[1] = max_xc[1] + p*dy[max_d];
		if(x[0]>=w || x[0]<0 || x[1]>=h || x[1]<0) break;
		if(*(mask+x[0]*h+x[1]) == 1) break;
		*(mask+x[0]*h+x[1]) = 1;
		p++;
	}

	p=1;
	while(1){
		x[0] = max_xc[0] - p*dx[max_d];
		x[1] = max_xc[1] - p*dy[max_d];
		if(x[0]>=w || x[0]<0 || x[1]>=h || x[1]<0) break;
		if(*(mask+x[0]*h+x[1]) == 1) break;
		*(mask+x[0]*h+x[1]) = 1;
		p++;
	}
			
}

