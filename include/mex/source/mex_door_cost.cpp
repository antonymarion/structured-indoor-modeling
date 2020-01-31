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
  
    double *ifsp0 = mxGetPr(prhs[0]); 
	double *ifsp1 = mxGetPr(prhs[1]);
	double *iwpf0 = mxGetPr(prhs[2]);
	double *iwpf1 = mxGetPr(prhs[3]);
    double *fsp0 = mxGetPr(prhs[4]);
	double *fsp1 = mxGetPr(prhs[5]);
	double *wp0 = mxGetPr(prhs[6]);
	double *wp1 = mxGetPr(prhs[7]);
	double lambda = (double)*mxGetPr(prhs[8]);
	int min_x = (int)*mxGetPr(prhs[9]);
	int max_x = (int)*mxGetPr(prhs[10]);
	int isDetail = (int)*mxGetPr(prhs[11]);

	double cf0, cf1, cw0, cw1;
	double cost;
	double min_cost = 1.0e6;
	double min_cost1, min_cost2;
	double min_x0, min_y0, min_x1, min_y1;
	

	int gNumY = (int)mxGetDimensions(prhs[4])[0];
	int gNumX = (int)mxGetDimensions(prhs[4])[1]; 
	//mexPrintf("%d %d\n",gNumY, gNumX);
    
    nlhs =4;
    plhs[0] = mxCreateDoubleMatrix(4, 1, mxREAL);
	double *X = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
	double *C = mxGetPr(plhs[1]);

	plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
	double *C1 = mxGetPr(plhs[2]);

	plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
	double *C2 = mxGetPr(plhs[3]);

	if(isDetail==0)
	{
		for(int x=min_x;x<=max_x;x++){
			for(int y=0;y<gNumY-1;y++){
				if(x<0 || x>=gNumX) continue;
				if(fsp0[x*gNumY+y]>0 && fsp1[x*gNumY+y]>0 && wp0[x*gNumY+y]==0 && wp1[x*gNumY+y]== 0){
				
					for(int xx=x+1;xx<=max_x;xx++){
						for(int yy=y+1;yy<gNumY;yy++){
							if(xx<0 || xx>=gNumX || yy <0 || yy>= gNumY) continue;
							if(fsp0[xx*gNumY+yy]>0 && fsp1[xx*gNumY+yy]>0 && wp0[xx*gNumY+yy]==0 && wp1[xx*gNumY+yy]== 0){			                                            
								int xx_ = xx+1;
								int yy_ = yy+1;
								if(xx_<0 || xx_>=gNumX+1 || yy_<0 || yy_>= gNumY+1) continue;
								cw0 = iwpf0[x*(gNumY+1)+y] + iwpf0[xx_*(gNumY+1)+yy_] - iwpf0[xx_*(gNumY+1)+y] - iwpf0[x*(gNumY+1)+yy_];
								cw1 = iwpf1[x*(gNumY+1)+y] + iwpf1[xx_*(gNumY+1)+yy_] - iwpf1[xx_*(gNumY+1)+y] - iwpf1[x*(gNumY+1)+yy_];
								cf0 = ifsp0[x*(gNumY+1)+y] + ifsp0[xx_*(gNumY+1)+yy_] - ifsp0[xx_*(gNumY+1)+y] - ifsp0[x*(gNumY+1)+yy_];
								cf1 = ifsp1[x*(gNumY+1)+y] + ifsp1[xx_*(gNumY+1)+yy_] - ifsp1[xx_*(gNumY+1)+y] - ifsp1[x*(gNumY+1)+yy_];
								cost = cw0+cw1-lambda*(cf0 + cf1);
								if(cost < min_cost){							
									min_x0 = x;
									min_y0 = y;
									min_x1 = xx;
									min_y1 = yy;
									min_cost = cost;
									min_cost1 = (cw0 - lambda*cf0)/((xx-x)*(yy-y));
									min_cost2 = (cw1 - lambda*cf1)/((xx-x)*(yy-y));
								}
							}
						}
					}
				
				}
			}
		}	
	}
	else
	{
		for(int x=min_x;x<=max_x;x++){
		for(int y=0;y<gNumY-1;y++){
				if(x<0 || x>=gNumX) continue;
				if(wp0[x*gNumY+y]==0 && wp1[x*gNumY+y]==0){
				
				for(int xx=x+1;xx<=max_x;xx++){
					for(int yy=y+1;yy<gNumY;yy++){
						if(xx<0 || xx>=gNumX || yy <0 || yy>= gNumY) continue;
						if(wp0[xx*gNumY+yy]== 0 && wp1[xx*gNumY+yy]== 0){			
							int xx_ = xx+1;
							int yy_ = yy+1;
							if(xx_<0 || xx_>=gNumX+1 || yy_<0 || yy_>= gNumY+1) continue;
							cw0 = iwpf0[x*(gNumY+1)+y] + iwpf0[xx_*(gNumY+1)+yy_] - iwpf0[xx_*(gNumY+1)+y] - iwpf0[x*(gNumY+1)+yy_];
							cw1 = iwpf1[x*(gNumY+1)+y] + iwpf1[xx_*(gNumY+1)+yy_] - iwpf1[xx_*(gNumY+1)+y] - iwpf1[x*(gNumY+1)+yy_];
							cf0 = ifsp0[x*(gNumY+1)+y] + ifsp0[xx_*(gNumY+1)+yy_] - ifsp0[xx_*(gNumY+1)+y] - ifsp0[x*(gNumY+1)+yy_];
							cf1 = ifsp1[x*(gNumY+1)+y] + ifsp1[xx_*(gNumY+1)+yy_] - ifsp1[xx_*(gNumY+1)+y] - ifsp1[x*(gNumY+1)+yy_];
							cost = cw0+cw1-lambda*(cf0 + cf1) - 100*(yy-y)*(xx-x)/double(gNumX*gNumY);
							if(cost < min_cost){
							
								min_x0 = x;
								min_y0 = y;
								min_x1 = xx;
								min_y1 = yy;
								min_cost = cost;
							}
					   }
					}
				}
				
			}
		}
		}	


	}
	*C = min_cost;
	*C1 = min_cost1;
	*C2 = min_cost2;
	*(X+0) = min_x0;
	*(X+1) = min_y0;
	*(X+2) = min_x1;
	*(X+3) = min_y1;
	
}
