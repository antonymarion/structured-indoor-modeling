// solve shortest path algorithm using Dijkstra method

#include <mex.h>
#include <math.h>
#include <iostream>
#include <queue>
#include <vector>
#include <functional>
using namespace std;

const double PI = 3.1415926;
struct State{
	int x;
	int y;

};

typedef pair<double, State> Node;

bool operator<(const State &a, const State &b) 
{ 
	if(a.x > b.x) return true;
	if(a.x == b.x && a.y > b.y) return true;
	return false;
} 

bool operator>(const State &a, const State &b) 
{ 
	return b < a; 
} 

bool operator==(const State &a, const State &b) 
{ 
	return a.x == b.x && a.y==b.y; 
} 

int sign(double x){
	
	if(x < 0) return -1;
	if(x == 0) return 0;
	return 1;

}
double min(double x, double y){
	if (x < y) return x;
	return y;
}

double norm(double v[2])
{
	return sqrt(v[0]*v[0]+v[1]*v[1]);
}

void cross(double out[3], const double v1[3], const double v2[3])
{
	out[0] = v1[1]*v2[2] - v1[2]*v2[1];
	out[1] = v1[2]*v2[0] - v1[0]*v2[2];
	out[2] = v1[0]*v2[1] - v1[1]*v2[0];
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	const double *cMap = mxGetPr(prhs[0]); // cost map
	const double *mask = mxGetPr(prhs[1]); // mask
	const double *fMap = mxGetPr(prhs[2]); // freespace evidence

	const double *start =mxGetPr(prhs[3]);
	const double *end = mxGetPr(prhs[4]);
	const double thresh = *mxGetPr(prhs[5]);
	const double alpha = *mxGetPr(prhs[6]);
	const int margin = (int)*mxGetPr(prhs[7]);
	const double sigma = *mxGetPr(prhs[8]);
	const double pmin = *mxGetPr(prhs[9]);
	const double fweight = *mxGetPr(prhs[10]);
	const double *nx = mxGetPr(prhs[11]);
	const double *ny = mxGetPr(prhs[12]);

	vector<double> kernel;
	double sum_kernel;

	const int h = mxGetDimensions(prhs[0])[0];
    const int w = mxGetDimensions(prhs[0])[1];
	


	vector<double> maskcopy(h*w);
	for(int i=0; i<h*w; i++) maskcopy[i] = mask[i];

	const int dx[] = {1, 0, -1, 0};
	const int dy[] = {0, 1, 0, -1};

	const int inout_index[] = {3, 0, 1, 2};


	double fincost;

	vector<double> closed(h*w);
	for (int i=1;i<h*w;i++) closed[i] = -1;
	vector<bool> done(h*w);
	for (int i=1;i<h*w;i++) done[i] = false;
	vector<int> edges_from_x(h*w);
	vector<int> edges_from_y(h*w);

	priority_queue <Node, vector<Node>, std::greater<Node> > open;
	State s = {(int)start[0], (int)start[1]};
	open.push(make_pair(0, s));
	closed[s.x*h + s.y] = 0;
	edges_from_x[s.x*h + s.y] = s.x;
	edges_from_y[s.x*h + s.y] = s.y;

	double dir_histgram[8];
	int count = 0;	
	bool flag2 = true;
	while(!open.empty()) {

		//if(count>4000) break;
		// pop out the node with minimum cost
		double cost = open.top().first;
		State cur = open.top().second;
		open.pop();
		done[cur.x*h+cur.y] = true;
		
		int pre_x = edges_from_x[cur.x*h+cur.y];
		int pre_y = edges_from_y[cur.x*h+cur.y];

		// for computationally efficiency
		int diffx = cur.x - pre_x;
		int diffy = cur.y - pre_y;
		int sdx = sign(diffx);
		int sdy = sign(diffy);

		if(sdx && !sdy){
			for(int d=1; d<sdx*diffx; d++) {maskcopy[(pre_x+d*sdx)*h + pre_y] = 1;}//
		}

		if(sdy && !sdx){
			for(int d=1; d<sdy*diffy; d++) maskcopy[pre_x*h + pre_y+d*sdy] = 1;
		}


		// judgement of the completion
		if (cur.x == (int)end[0] && cur.y == (int)end[1]) 
		{
			fincost = cost;
			break;
		}

		// Path only for Manhattan World Directions
		for (int d=0; d<4; d++) {

			// for computationally efficiency
			if(cur.x - pre_x != 0 && cur.y - pre_y == 0 && (d == 0 || d == 2)) continue; // restrict edges in the vertical direction
			if(cur.y - pre_y != 0 && cur.x - pre_x == 0 && (d == 1 || d == 3)) continue; // restrict edges in the horizontal direction
			if(cur.x == start[0] && cur.y == start[1] && (end[0]-start[0])*(end[0]-start[0]) > 0 && (d==1||d==3)) continue;
			if(cur.x == start[0] && cur.y == start[1] && (end[1]-start[1])*(end[1]-start[1]) > 0 && (d==0||d==2)) continue;

			int p = 1;
			double costOn = 0;
			double costPath = alpha;
			double sum_n[2] = {};

			while (1) {
				int new_x = cur.x + p*dx[d];
				int new_y = cur.y + p*dy[d];

				double new_nx=0;
				double new_ny=0;

				p++;
				bool flag = true;
				if (new_x >= 0 && new_x < w && new_y >= 0 && new_y < h)
				{		
					double cPoint = 0;
					double cFree = 0;

					cPoint = cMap[new_x*h+new_y]; 
					cFree = fMap[new_x*h+new_y];
					new_nx = nx[new_x*h+new_y];
					new_ny = ny[new_x*h+new_y];

					double nn[2];
					nn[0] = new_nx;
					nn[1] = new_ny;

					double nm = norm(nn);
					double dir[3] = {dx[d], dy[d], 0};
					double cur_cross[3];
					
					if(nm>0)
					{
					double cur_n[3] = {nn[0]/nm, nn[1]/nm, 0};
					cross(cur_cross, dir, cur_n);
					}
					else
					{
					cur_cross[2] = 1; 
					}
												

				
					double weight = (1-cPoint + fweight*cFree)/(1+fweight);
					if(cPoint == 0){
						weight = weight*4;
					}

					costOn = costOn + weight;

					if(mask[new_x*h+new_y] == 1)
					{
						costOn += 1.0e6;
						break;
					}


					double penalty = 0;
					
					if(p-1 < pmin && (cur.x != start[0] || cur.y != start[1]) && (new_x != end[0] || new_y != end[1]))
					{
						penalty += 1.0e6;
					}

					//if(cur_cross[2] < -0.90)
					//{
					//	costOn += 1.0e6;
					//	break;
					//}
					

					double new_cost = cost + costOn + penalty;
					if((new_x != end[0] || new_y != end[1]) & (cur.x != start[0] || cur.y != start[1])) new_cost += costPath;
					State new_state = {new_x, new_y};
					count = count + 1;
					if((closed[new_state.x*h+new_state.y] == -1 || new_cost < closed[new_state.x*h+new_state.y]) && done[new_state.x*h+new_state.y] == false){
						closed[new_state.x*h+new_state.y] = new_cost;
						edges_from_x[new_state.x*h+new_state.y] = cur.x;
						edges_from_y[new_state.x*h+new_state.y] = cur.y;
						open.push(make_pair(new_cost, new_state));
					}
				}
				else
				{
					break;
				}
			}
		}
	}




	// inversely trace the shortest path
	State e = {(int)end[0], (int)end[1]};
	vector <int> index;
	
	double sum_path_length = 0;
	while(1){
		index.push_back(e.x*h+e.y);
		State next;
		next.x = edges_from_x[e.x*h+e.y];
		next.y = edges_from_y[e.x*h+e.y];
		sum_path_length += sqrt(double((next.x-e.x)*(next.x-e.x) + (next.y-e.y)*(next.y-e.y)));
		e = next;
		if(next.x == start[0] && next.y == start[1]) {
			index.push_back(next.x*h+next.y);
			break;
		}

	}

	nlhs = 1;
    plhs[0] = mxCreateDoubleMatrix(index.size(), 1, mxREAL);
	double *shortest_path = mxGetPr(plhs[0]);
	for (int i=0; i<index.size();i++)
	*(shortest_path+i) = index[i];

	plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
	double *min_cost = mxGetPr(plhs[1]);
	*min_cost = fincost;

	plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
	double *length_path = mxGetPr(plhs[2]);
	*length_path = sum_path_length;

	plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
	double *cost_count = mxGetPr(plhs[3]);
	*cost_count = count;
	
}

