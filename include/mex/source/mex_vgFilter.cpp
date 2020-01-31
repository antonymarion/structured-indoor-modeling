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

	vector<vector<unsigned long long> > index_list;
	/*vector<unsigned long long> list_to_voxel;*/
	vector<unsigned int> voxel_to_list(gNum[0]*gNum[1]*gNum[2]);
	double x[3];
	unsigned int idx, idy, idz;
	unsigned long long list_id;
	for(int i=0; i<m; i++)
	{
		x[0] = P[3*i+0];
		x[1] = P[3*i+1];
		x[2] = P[3*i+2];
		idx = floor(gNum[0]*(x[0]-b0[0])/bSize[0]);
		idy = floor(gNum[1]*(x[1]-b0[1])/bSize[1]);
		idz = floor(gNum[2]*(x[2]-b0[2])/bSize[2]);
		unsigned long long index = idz*gNum[0]*gNum[1]+idx*gNum[1]+idy;
		if(idx>=0 && idx<gNum[0] && idy>=0 && idy<gNum[1] && idz>=0 && idz<gNum[2])
		{
			if(voxel_to_list[index] == 0){
				vector<unsigned long long> sub_index_list;
				sub_index_list.push_back(i);
				index_list.push_back(sub_index_list);
				voxel_to_list[index] = index_list.size();
				//list_to_voxel.push_back(index);
			}
			else
			{
				list_id = voxel_to_list[index];
				index_list[list_id-1].push_back(i);
			}
		}
	}

	unsigned int num_point;
	unsigned int point_index, medoid_index;
	double d, medoid_distance;
	nlhs = 1;
    plhs[0] = mxCreateDoubleMatrix(index_list.size(), 1, mxREAL);
	double *filtered_index = mxGetPr(plhs[0]);

	for(int i=0; i<index_list.size(); i++)
	{
		num_point = index_list[i].size();
		vector<double> xlist(num_point);
		vector<double> ylist(num_point);
		vector<double> zlist(num_point);
		
		// strage point
		for(int j=0; j<num_point; j++)
		{
			point_index = index_list[i][j];
			xlist[j] = P[3*point_index+0];
			ylist[j] = P[3*point_index+1];
			zlist[j] = P[3*point_index+2];
		}

		// compute medoids
		medoid_index = -1;
		medoid_distance = DBL_MAX;
		for(int j=0; j<num_point; j++)
		{
			d = 0;
			for(int k=0; k<num_point; k++)
			{
				if(k!=j) d += (xlist[j]-xlist[k])*(xlist[j]-xlist[k]) + (ylist[j]-ylist[k])*(ylist[j]-ylist[k]) + (zlist[j]-zlist[k])*(zlist[j]-zlist[k]);
			}

			if(d < medoid_distance)
			{
				medoid_index = j;
				medoid_distance = d;
			}
		}

		 *(filtered_index+i) = index_list[i][medoid_index]+1;
	}
    
    
	


	
	
}
