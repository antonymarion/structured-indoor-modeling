#include <mex.h>
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	const int width = (int)*mxGetPr(prhs[0]);
	const int height = (int)*mxGetPr(prhs[1]);
	double *mask = mxGetPr(prhs[2]);    
    
	int index, flag;
    nlhs = 1;
    plhs[0] = mxCreateDoubleMatrix(width*height, 1, mxREAL);
	double *boundary = mxGetPr(plhs[0]);
	int x0, x1, y0, y1;

	// boundary 
	/* x=0, x=width-1*/
	for(int y=0; y < height; ++y)
	{
		x0 = 0;
		index = x0 * height + y;
		flag = (int)mask[index];
		if (flag==1)
		{	
			*(boundary+index) = 1;	
		}

		x0 = width-1;
		index = x0 * height + y;
		flag = (int)mask[index];
		if (flag==1)
		{	
			*(boundary+index) = 1;	
		}
	}


	/* y=0, y=height-1*/
	for(int x=0; x < width; ++x)
	{
		y0 = 0;
		index = x * height + y0;
		flag = (int)mask[index];
		if (flag==1)
		{	
			*(boundary+index) = 1;	
		}

		y0 = height-1;
		index = x * height + y0;
		flag = (int)mask[index];
		if (flag==1)
		{	
			*(boundary+index) = 1;	
		}
	}


	for (int y = 1; y < height - 1; ++y) {
		for (int x = 1; x < width - 1; ++x) {
			index = x * height + y;
			flag = (int)mask[index];
			if (flag==0)
			{				
				continue;
			}
			
			if ((!mask[index - 1]) ||
			(!mask[index + 1]) ||
			(!mask[index - height]) ||
			(!mask[index + height])) {
			*(boundary+index) = 1;
			}
				
		}
	}	
}

