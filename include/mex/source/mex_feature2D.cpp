#include <mex.h>
#include <math.h>

int sign(double x)
{
	if(x >= 0) 
		return 1;
	else
		return -1;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  
    double *gNum = mxGetPr(prhs[0]); // grid of the bounding box
	bool *FS = mxGetLogicals(prhs[1]); // FS evidence  (size: gNum(1)*gNum(2))
    double *V = mxGetPr(prhs[2]); // directions of the ray (in float space) (size: 3*num_vec)
	int gNumX = (int)gNum[0];
	int gNumY = (int)gNum[1];
	int num_vec = mxGetDimensions(prhs[2])[1];
    
    nlhs = 1;
    plhs[0] = mxCreateDoubleMatrix(gNumX*gNumY*num_vec, 1, mxREAL);
	double *F = mxGetPr(plhs[0]);
	double kX[2];
	double v[2];
	int step[2];
	int X[2];
	double tDelta[2];
	double tMax[2];
	
	for (int x=0; x<gNumX; x++)
	{
		for (int y=0; y<gNumY; y++)
		{
			if (FS[y*gNumX+x] == 1)
			{
				for(int i=0; i<num_vec; i++)
				{

					v[0] = V[2*i+0];
					v[1] = V[2*i+1];

					step[0] = sign(v[0]);
					step[1] = sign(v[1]);

					X[0] = x;
					X[1] = y;


					// how far away along the ray we must move
					if(v[0]*v[0] < 1.0e-12){
						tDelta[0] = 1.0e12;
						tMax[0] = 1.0e12;
					}
					else{
						tDelta[0] = abs(1/v[0]);
						tMax[0] = abs(1/v[0]);
					}

					if(v[1]*v[1] < 1.0e-12){
						tDelta[1] = 1.0e12;
						tMax[1] = 1.0e12;
					}
					else{
						tDelta[1] = abs(1/v[1]);
						tMax[1] = abs(1/v[1]);
					}


					int count = 0;
					while (1)
					{			
						if(X[1]*gNumX+X[0] < 0 || X[1]*gNumX+X[0] >= gNumX*gNumY)
						{
							F[i*gNumX*gNumY + (y*gNumX+x)] = sqrt(double((X[0]-x)*(X[0]-x)+(X[1]-y)*(X[1]-y)));
							break;
						}
			
						if (FS[X[1]*gNumX+X[0]]==0)
						{
							F[i*gNumX*gNumY + (y*gNumX+x)] = sqrt(double((X[0]-x)*(X[0]-x)+(X[1]-y)*(X[1]-y)));
							break;
						}

						if (tMax[0] < tMax[1])
						{
								tMax[0] = tMax[0] + tDelta[0];
								X[0] = X[0] + step[0];		
								
						}
						else
						{
								tMax[1] = tMax[1] + tDelta[1];
								X[1] = X[1] + step[1];	
							
						}
						count = count + 1;
					}



				}
			}
		}
	}

	
		
	
}
