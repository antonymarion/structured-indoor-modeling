#include <mex.h>
#include <math.h>

double min(double x, double y)
{
	if(x > y)
		return y;
	else
		return x;
	
}
bool IsVisible(const int width, const int height, const int x0, const int y0, const int x1, const int y1, const double *mask)
{
	double v[2];
	int x, y;
	double step_x, step_y;
	v[0] = x1 - x0;
	v[1] = y1 - y0;
	double t = min(width, height);
	const int num_steps = min(t, 2*floor(sqrt(v[0]*v[0]+v[1]*v[1])));
	step_x = v[0]/num_steps;
	step_y = v[1]/num_steps;

	for (int i = 1; i < num_steps; i++)
	{
		x = (int) floor(x0 + i * step_x);
		y = (int) floor(y0 + i * step_y);
		const int index = x*height + y;
		if(index < 0) return false;
		if(index >= width*height) return false;
		if (!mask[index])
		{
			return false;
		}
	}

	return true;	
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    const int width = (int)*mxGetPr(prhs[0]);
    const int height = (int)*mxGetPr(prhs[1]);
    double *mask = mxGetPr(prhs[2]);
    double *seed = mxGetPr(prhs[3]);
    double *target = mxGetPr(prhs[4]); 
	int num_seed = mxGetDimensions(prhs[3])[0];
    int num_target = mxGetDimensions(prhs[4])[0];


    nlhs = 1;
	plhs[0] = mxCreateLogicalMatrix(num_target*num_seed, 1);
    bool *feature =  mxGetLogicals(plhs[0]);
	bool visibility;
	int x0, y0, x1, y1;
	
    for (int j = 0; j < num_target; j++)
	{
		x0 = (int)target[j]-1;
		y0 = (int)target[j+num_target]-1;
		for (int i = 0; i < num_seed; i++)
		{
			x1 = (int)seed[i]-1;
			y1 = (int)seed[num_seed+i]-1;
			visibility = IsVisible(width, height, x0, y0, x1, y1, mask);
			*(feature + num_seed*j + i) = visibility;
		}
	}





}
