#include <mex.h>
#include <math.h>
#include <deque>
#include <vector>
using namespace std;

typedef pair<int, int> Node;

int sign(double x)
{
	if(x >= 0) 
		return 1;
	else
		return -1;
}

int min(int a, int b){
	if(a >= b)
		return b;
	else
		return a;
}

int max(int a, int b){
	if(a >= b)
		return a;
	else
		return b;
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
	vector<double> used_list(gNumX*gNumY);
	double *maskOut = mxGetPr(plhs[0]);

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
		vx = vx/nv;
		vy = vy/nv;
		double count = 0;
		int num_step = 1000;
		double step = nv/num_step;
		
		for(int j=0;j<=num_step;j++){
			x = floor(x0 + j*step*vx);
			y = floor(y0 + j*step*vy);
			mask[x*gNumY+y] = 1;
		}
	}

	int xstart, ystart;
	int yset0, yset1;
	// find starting point
	int count;
	for(int i=0;i<pNum-1;i++){
		x0 = xlist[i]-1;
		x1 = xlist[i+1]-1;
		int minx = min(x0,x1);
		int maxx = max(x0,x1);
		for(int x=minx; x<=maxx;x++)
		{
			count = 0;
			for(int y=0;y<gNumY;y++){
			
				if(mask[x*gNumY+y] == 1 && count == 1){
					count = 0;
					if(abs(yset0-y)>1){
					yset1 = y;
					//mexPrintf("x1 y1 %d %d\n", x, yset1);
					xstart = x+1;
					ystart = (yset0+yset1)/2;
					count = 2;
					}				
				}

				if(mask[x*gNumY+y] == 1 && count == 0){
					yset0 = y;
					//mexPrintf("x0 y0 %d %d\n", x, yset0);
					count = 1;
				}		

				if(count == 2) break;
			}
			if(count == 2) break;
		}
		if(count == 2) break;
	}
	//for(int i=0;i<gNumX*gNumY;i++) *(maskOut+i) = mask[i];
	//mexPrintf("%d %d %d\n", count, xstart, ystart);
	int test_count = 0;
	if(count < 2){for(int i=0;i<gNumX*gNumY;i++) *(maskOut+i) = mask[i];}
	else{
		deque<Node> nodelist;
		nodelist.push_back(make_pair(xstart-1,ystart));
		while(!nodelist.empty())
		{
			test_count++;
			if(test_count > 1000) break;
			int x = nodelist.front().first;
			int y = nodelist.front().second;
			//mexPrintf("%d %d A: %d %d\n", test_count, nodelist.size(), x, y);
			nodelist.pop_front();
			
			//int xx = nodelist.front().first;
			//int yy = nodelist.front().second;
			/*mexPrintf("%d %d B: %d %d\n", test_count, nodelist.size(), x, y);*/
			mask2[x*gNumY+y] = 1;
			//left
			int cnt = 0;
			int xleft, xright;
			xleft = x;
			xright = x;
			while(1){
			
				cnt++;
				int xx = x - cnt; 
				int yy = y;
				//mexPrintf("%d %d\n", xx, yy);
				if(mask[xx*gNumY+y] == 1) {xleft = x - cnt + 1; /*mexPrintf("LEFT %d %d\n", xx, yy);*/break;}
				mask2[xx*gNumY+y] = 1;
				
			}

			cnt = 0;
			while(1){
				cnt++;
				int xx = x + cnt; 
				int yy = y;
				//mexPrintf("%d %d\n", xx, yy);
				if(mask[xx*gNumY+y] == 1) {xright = x + cnt - 1; /*mexPrintf("RIGHT %d %d\n", xx, yy);*/break;}
				mask2[xx*gNumY+y] = 1;

			}

			// scan from left to right (y-1)
			for(int p = xleft; p <= xright; p++){
				if(p != xright && mask[p*gNumY+y-1] == 0 && mask[(p+1)*gNumY+y-1] == 1 && used_list[p*gNumY+y-1] == 0)
				{
				nodelist.push_back(make_pair(p,y-1));
				used_list[p*gNumY+y-1] = 1;
				//mexPrintf("PUT: %d %d\n", p, y-1);
				}

				if(p == xright && mask[p*gNumY+y-1] == 0 && used_list[xright*gNumY+y-1] == 0) 
				{
					nodelist.push_back(make_pair(xright,y-1));
					used_list[xright*gNumY+y-1] = 1;
				}

				if(p != xright, mask[p*gNumY+y+1] == 0 && mask[(p+1)*gNumY+y+1] == 1 && used_list[p*gNumY+y+1] == 0)
				{
					nodelist.push_back(make_pair(p,y+1));
					used_list[p*gNumY+y+1] = 1;
				//mexPrintf("PUT: %d %d\n", p, y+1);
				}

				if(p == xright && mask[p*gNumY+y+1] == 0 && used_list[xright*gNumY+y+1] == 0)
				{
					nodelist.push_back(make_pair(xright,y+1));
					used_list[xright*gNumY+y+1] = 1;
				}
			}
		}
		for(int i=0;i<gNumX*gNumY;i++){ *(maskOut+i) = mask2[i] + mask[i];};
	}
	
}
