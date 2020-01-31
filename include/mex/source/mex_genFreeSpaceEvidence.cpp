#include <mex.h>
#include <stdlib.h>
#include <cmath>
#include <stdio.h>
#include <time.h>
#include <iostream>
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
    double *P = mxGetPr(prhs[0]); 
	double *C = mxGetPr(prhs[1]);
	double *CAMLIST = mxGetPr(prhs[2]);
	double *gNum = mxGetPr(prhs[3]);
	double *leaf_size = mxGetPr(prhs[4]);
	double *b0 = mxGetPr(prhs[5]);
	double ratio = *mxGetPr(prhs[6]);
	int m = (int)mxGetDimensions(prhs[0])[1];
	double x[3];
	unsigned int idx, idy, idz;

    nlhs = 1;
	plhs[0] = mxCreateDoubleMatrix(gNum[0]*gNum[1]*gNum[2], 1, mxREAL);
    double *freespace_evidence =  mxGetPr(plhs[0]);

	// initialization
	int X0[3];
	int X[3];		
	double p0[3], p1[3];
	int Y[3];
	int cid;
	double v[3];
	int step[3];
	double tDelta[3];
	double tMax[3];

	for (int i=0;i<m;i++)
	{
		p1[0] = P[3*i+0];
		p1[1] = P[3*i+1];
		p1[2] = P[3*i+2];

		cid = C[i]-1;
		
		p0[0] = CAMLIST[3*cid+0];
		p0[1] = CAMLIST[3*cid+1];
		p0[2] = CAMLIST[3*cid+2];
		
		X[0] = floor((p0[0]-b0[0])/leaf_size[0]);
		X[1] = floor((p0[1]-b0[1])/leaf_size[1]);
		X[2] = floor((p0[2]-b0[2])/leaf_size[2]);

		Y[0] = floor((p1[0]-b0[0])/leaf_size[0]);
		Y[1] = floor((p1[1]-b0[1])/leaf_size[1]);                               
		Y[2] = floor((p1[2]-b0[2])/leaf_size[2]);

		v[0] = p1[0] - p0[0];
		v[1] = p1[1] - p0[1];
		v[2] = p1[2] - p0[2];

		step[0] = sign(v[0]);
		step[1] = sign(v[1]);
		step[2] = sign(v[2]);
			
		tDelta[0] = abs(leaf_size[0]/v[0]);
		tDelta[1] = abs(leaf_size[1]/v[1]);
		tDelta[2] = abs(leaf_size[2]/v[2]);

		tMax[0] = abs((0.5*(1-step[0])*leaf_size[0]-(b0[0] + leaf_size[0]*(X[0]+1) - p0[0]))/v[0]);
		tMax[1] = abs((0.5*(1-step[1])*leaf_size[1]-(b0[1] + leaf_size[1]*(X[1]+1) - p0[1]))/v[1]);
		tMax[2] = abs((0.5*(1-step[2])*leaf_size[2]-(b0[2] + leaf_size[2]*(X[2]+1) - p0[2]))/v[2]);
		
	    
		int count = 0;
		unsigned long long index;
		vector<int> xlist;
		vector<int> ylist;
		vector<int> zlist;
		
		while (1)
		{		


			if (X[0] == Y[0] && X[1] == Y[1] && X[2] == Y[2]) break;

			if (X[0] < 0 || X[0] >= gNum[0] || X[1] < 0 || X[1] >= gNum[1] || X[2] < 0 || X[2] >= gNum[2]) break;

			if (tMax[0] < tMax[1])
			{
				if (tMax[0] < tMax[2])
				{
					tMax[0] = tMax[0] + tDelta[0];
					X[0] = X[0] + step[0];
					index = X[2]*gNum[0]*gNum[1]+X[0]*gNum[1]+X[1];
					if(X[0]>=0 && X[0]<gNum[0] && X[1]>=0 && X[1]<gNum[1] && X[2]>=0 && X[2]<gNum[2])
					{
					xlist.push_back(X[0]);
					ylist.push_back(X[1]);
					zlist.push_back(X[2]);
					}
					
				}
				else
				{
					tMax[2] = tMax[2] + tDelta[2];
					X[2] = X[2] + step[2];
					index = X[2]*gNum[0]*gNum[1]+X[0]*gNum[1]+X[1];
					if(X[0]>=0 && X[0]<gNum[0] && X[1]>=0 && X[1]<gNum[1] && X[2]>=0 && X[2]<gNum[2])
					{
					xlist.push_back(X[0]);
					ylist.push_back(X[1]);
					zlist.push_back(X[2]);
					}
					
				}
			}
			else
			{
				if (tMax[1] < tMax[2])
				{
					tMax[1] = tMax[1] + tDelta[1];
					X[1] = X[1] + step[1];
					index = X[2]*gNum[0]*gNum[1]+X[0]*gNum[1]+X[1];
					if(X[0]>=0 && X[0]<gNum[0] && X[1]>=0 && X[1]<gNum[1] && X[2]>=0 && X[2]<gNum[2])
					{
					xlist.push_back(X[0]);
					ylist.push_back(X[1]);
					zlist.push_back(X[2]);
					}
				}
				else
				{
					tMax[2] = tMax[2] + tDelta[2];
					X[2] = X[2] + step[2];
					index = X[2]*gNum[0]*gNum[1]+X[0]*gNum[1]+X[1];
					if(X[0]>=0 && X[0]<gNum[0] && X[1]>=0 && X[1]<gNum[1] && X[2]>=0 && X[2]<gNum[2])
					{
					xlist.push_back(X[0]);
					ylist.push_back(X[1]);
					zlist.push_back(X[2]);					
					}
				}

			}
		}

		int size_list = xlist.size();

		for(int j=0;j<ratio*xlist.size();j++)
		{

			index = zlist[j]*gNum[0]*gNum[1]+xlist[j]*gNum[1]+ylist[j];
			*(freespace_evidence+index) = *(freespace_evidence+index) + 1;

		}
		
	}   
}
