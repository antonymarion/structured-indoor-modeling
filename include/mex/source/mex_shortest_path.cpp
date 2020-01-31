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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	const double *sList = mxGetPr(prhs[0]); // startlist
	const double *target_mask = mxGetPr(prhs[1]); // mask
	const double *free_mask = mxGetPr(prhs[2]);
	const int h = mxGetDimensions(prhs[1])[0];
    const int w = mxGetDimensions(prhs[1])[1];
	const int numseed = mxGetDimensions(prhs[0])[1];

	const int dx[] = {1, 0, -1, 0};
	const int dy[] = {0, 1, 0, -1};


	int start[2];
	int end[2];
	int mincost = INT_MAX;
	int min_end[2];
	int min_start[2];
	vector<int> min_edges_from_x(h*w);
	vector<int> min_edges_from_y(h*w);


	vector<bool> finished(h*w);
	for(int j=0; j<numseed; j++)
	{
		
		//mexPrintf("%d\n",j);
		start[0] = sList[2*j+0];
		start[1] = sList[2*j+1];
		int fincost;
		vector<double> closed(h*w);
		for (int i=1;i<h*w;i++) closed[i] = -1;
		vector<bool> done(h*w);
		for (int i=1;i<h*w;i++) done[i] = false;
		vector<int> edges_from_x(h*w);
		vector<int> edges_from_y(h*w);
	
		priority_queue <Node, vector<Node>, greater<Node> > open;
		State s = {(int)start[0], (int)start[1]};
		open.push(make_pair(0, s));
		closed[s.x*h + s.y] = 0;

		//mexPrintf("%d-th start point\n", j);
		int count = 0;
		while(!open.empty()) {
			// pop out the node with minimum cost
			double new_cost = open.top().first + 1;
			State cur = open.top().second;
			open.pop();
			done[cur.x*h+cur.y] = true;

			// judgement of the completion
			if (target_mask[cur.x*h + cur.y] == 1) 
			{
				fincost = new_cost;
				end[0] = cur.x;
				end[1] = cur.y;
				break;
			}

			// Path only for Manhattan World Directions
			for (int d=0; d<4; d++) {

				int new_x = cur.x + dx[d];
				int new_y = cur.y + dy[d];
		
				bool flag = true;
				if (new_x >= 0 && new_x < w && new_y >= 0 && new_y < h)
				{
					if(free_mask[new_x*h + new_y] == 1){
						State new_state = {new_x, new_y};
						if((closed[new_state.x*h+new_state.y] == -1 || new_cost < closed[new_state.x*h+new_state.y]) && done[new_state.x*h+new_state.y] == false && new_cost < mincost){
							count++;
							closed[new_state.x*h+new_state.y] = new_cost;
							edges_from_x[new_state.x*h+new_state.y] = cur.x;
							edges_from_y[new_state.x*h+new_state.y] = cur.y;
							open.push(make_pair(new_cost, new_state));
						}
					}
				}			
			}
		}

		if(!open.empty())
		{
			mincost = fincost;
			min_edges_from_x = edges_from_x;
			min_edges_from_y = edges_from_y;
			min_end[0] = end[0];
			min_end[1] = end[1];
			min_start[0] = start[0];
			min_start[1] = start[1];
		}
	}


	// inversely trace the shortest path
	State e = {(int)min_end[0], (int)min_end[1]};
	/*vector <int> index;*/
	

	while(1){
		finished[e.x*h+e.y] = true;
		State next;
		next.x = min_edges_from_x[e.x*h+e.y];
		next.y = min_edges_from_y[e.x*h+e.y];
		e = next;
		if(next.x == min_start[0] && next.y == min_start[1]) {
			finished[e.x*h+e.y] = true;
			break;
		}
	}
	
	

	/*nlhs = 1;
    plhs[0] = mxCreateDoubleMatrix(index.size(), 1, mxREAL);
	double *shortest_path = mxGetPr(plhs[0]);
	for (int i=0; i<index.size();i++)
	*(shortest_path+i) = index[i];*/
	
	nlhs = 1;
    plhs[0] = mxCreateDoubleMatrix(h*w, 1, mxREAL);
	double *shortest_path = mxGetPr(plhs[0]);
	for (int i=0; i<h*w;i++)
	*(shortest_path+i) = finished[i];
}

