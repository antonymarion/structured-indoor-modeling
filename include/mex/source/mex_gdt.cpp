/* c++ mex implementation of geodesic distance transform for joint geodesic distance upsampling*/

#include "mex.h"
#include <iostream>
#include <map>
#include <cmath>    
using namespace std;

// return pointer index of matrix[s][5]
double access(const double *array, int s, int k, int h, int w)
{
    int id = k*h*w + s;
    return *(array+id);
}

std::map<int, double> minIDandVal(double *val, int m)
{
    std::map<int, double> r;

    double max_val = val[0];
    int max_id = 0;

    for(int i = 1; i < m; i++)
    {
        if(val[i] < max_val)
        {
            max_val = val[i];
            max_id = i;
        }
    }

    r.insert(std::map<int, double>::value_type(max_id, max_val));
    return r;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *M = mxGetPr(prhs[0]);
    double *L = mxGetPr(prhs[1]);
    double *dIf = mxGetPr(prhs[2]);
    double *dIb = mxGetPr(prhs[3]);
    double r = *mxGetPr(prhs[4]);
    double lambda = *mxGetPr(prhs[5]);
    int maxIter = (int)*mxGetPr(prhs[6]);



    int i, j, h, w;
    int s0, s1, s2, s3, s4;
    double d0, d1, d2, d3, d4;
    double dI[5];
    h = mxGetDimensions(prhs[0])[0];
    w = mxGetDimensions(prhs[0])[1];


    // L2 distance from neighbor pixels
    double dx[5] = {0, sqrt(2.0), 1.0, sqrt(2.0), 1.0};
    double vfx[5] = {0, -1, 0, 1, -1};
    double vfy[5] = {0, -1, -1, -1, 0};
    double vbx[5] = {0, 1, 0, -1, 1};
    double vby[5] = {0, 1, 1, 1, 0};


    nlhs = 2;
    plhs[0] = mxCreateDoubleMatrix(h, w, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(h, w, mxREAL);


    double *Mout = mxGetPr(plhs[0]);
    double *Lout = mxGetPr(plhs[1]);



    for(int iter = 0; iter < maxIter; iter++)
    {



        // forward scan update
        // case j = 1;
        j = 0;
        for(i = 1; i < w; i++)
        {

            s0 = i*h + j;
            s4 = (i+vfx[4])*h + j + vfy[4];

            dI[0] = access(dIf, s0, 0, h, w);
            dI[4] = access(dIf, s4, 4, h, w);

            d0 = *(M+s0);
            d4 = 1/r*dx[4] + lambda*dI[4] + *(M+s4);

            if (d4 < d0)
            {
                *(M+s0) = d4;
                *(L+s0) = *(L+s4);
            }

        }


        for(j = 1; j < h; j++)
        {
            // case i = 1
            i = 0;
            s0 = i*h + j;
            s2 = (i+vfx[2])*h + j + vfy[2];
            s3 = (i+vfx[3])*h + j + vfy[3];

            dI[0] = access(dIf, s0, 0, h, w);
            dI[2] = access(dIf, s2, 2, h, w);
            dI[3] = access(dIf, s3, 3, h, w);

            d0 = *(M+s0);
            d2 = 1/r*dx[2] + lambda*dI[2] + *(M+s2);
            d3 = 1/r*dx[3] + lambda*dI[3] + *(M+s3);

            double val[4] = {d0, d2, d3};
            map<int, double> idval;
            idval = minIDandVal(val, 3);
            map<int, double>::iterator it = idval.begin();

            switch((*it).first)
            {
            case 0: break;
            case 1: *(M+s0) = d2; *(L+s0) = *(L+s2); break;
            case 2: *(M+s0) = d3; *(L+s0) = *(L+s3); break;
            }

            for(i = 1; i<w-1; i++)
            {
                s0 = i*h + j;
                s1 = (i+vfx[1])*h + j + vfy[1];
                s2 = (i+vfx[2])*h + j + vfy[2];
                s3 = (i+vfx[3])*h + j + vfy[3];
                s4 = (i+vfx[4])*h + j + vfy[4];

                dI[0] = access(dIf, s0, 0, h, w);
                dI[1] = access(dIf, s1, 1, h, w);
                dI[2] = access(dIf, s2, 2, h, w);
                dI[3] = access(dIf, s3, 3, h, w);
                dI[4] = access(dIf, s4, 4, h, w);


                d0 = *(M+s0);
                d1 = 1/r*dx[1] + lambda*dI[1] + *(M+s1);
                d2 = 1/r*dx[2] + lambda*dI[2] + *(M+s2);
                d3 = 1/r*dx[3] + lambda*dI[3] + *(M+s3);
                d4 = 1/r*dx[4] + lambda*dI[4] + *(M+s4);

                double val[5] = {d0, d1, d2, d3, d4};
                map<int, double> idval;
                idval = minIDandVal(val, 5);
                map<int, double>::iterator it = idval.begin();
                switch((*it).first)
                {
                case 0: break;
                case 1: *(M+s0) = d1; *(L+s0) = *(L+s1); break;
                case 2: *(M+s0) = d2; *(L+s0) = *(L+s2); break;
                case 3: *(M+s0) = d3; *(L+s0) = *(L+s3); break;
                case 4: *(M+s0) = d4; *(L+s0) = *(L+s4); break;
                }

            }

            // case i = w
            i = w-1;
            s0 = i*h + j;
            s1 = (i+vfx[1])*h + j + vfy[1];
            s2 = (i+vfx[2])*h + j + vfy[2];
            s4 = (i+vfx[4])*h + j + vfy[4];

            dI[0] = access(dIf, s0, 0, h, w);
            dI[1] = access(dIf, s1, 1, h, w);
            dI[2] = access(dIf, s2, 2, h, w);
            dI[4] = access(dIf, s4, 4, h, w);

            d0 = *(M+s0);
            d1 = 1/r*dx[1] + lambda*dI[1] + *(M+s1);
            d2 = 1/r*dx[2] + lambda*dI[2] + *(M+s2);
            d4 = 1/r*dx[4] + lambda*dI[4] + *(M+s4);

            val[0] = d0; val[1] =  d1; val[2] = d2; val[3] =  d4;
            idval = minIDandVal(val, 4);
            it = idval.begin();
            switch((*it).first)
            {
            case 0: break;
            case 1: *(M+s0) = d1; *(L+s0) = *(L+s1); break;
            case 2: *(M+s0) = d2; *(L+s0) = *(L+s2); break;
            case 3: *(M+s0) = d4; *(L+s0) = *(L+s4); break;
            }
        }



        // backward scan update
        // case j = h;
        j = h-1;

        for(i = w-2; i >= 0; i--)
        {

            s0 = i*h + j;
            s4 = (i+vbx[4])*h + j + vby[4];
            dI[0] = access(dIb, s0, 0, h, w);
            dI[4] = access(dIb, s4, 4, h, w);

            d0 = *(M+s0);
            d4 = 1/r*dx[4] + lambda*dI[4] + *(M+s4);
            if (d4 < d0)
            {
                *(M+s0) = d4;
                *(L+s0) = *(L+s4);
            }
        }

        for(j = h-2; j >= 0; j--)
        {
            // case i = w
            i = w-1;
            s0 = i*h + j;
            s2 = (i+vbx[2])*h + j + vby[2];
            s3 = (i+vbx[3])*h + j + vby[3];
            dI[0] = access(dIb, s0, 0, h, w);
            dI[2] = access(dIb, s2, 2, h, w);
            dI[3] = access(dIb, s3, 3, h, w);


            d0 = *(M+s0);
            d2 = 1/r*dx[2] + lambda*dI[2] + *(M+s2);
            d3 = 1/r*dx[3] + lambda*dI[3] + *(M+s3);

            double val[3] = {d0, d2, d3};
            map<int, double> idval;
            idval = minIDandVal(val, 3);
            map<int, double>::iterator it = idval.begin();
            switch((*it).first)
            {
            case 0: break;
            case 1: *(M+s0) = d2; *(L+s0) = *(L+s2); break;
            case 2: *(M+s0) = d3; *(L+s0) = *(L+s3); break;
            }

            for(i = w-2; i >= 1; i--)
            {
                s0 = i*h + j;
                s1 = (i+vbx[1])*h + j + vby[1];
                s2 = (i+vbx[2])*h + j + vby[2];
                s3 = (i+vbx[3])*h + j + vby[3];
                s4 = (i+vbx[4])*h + j + vby[4];

                dI[0] = access(dIb, s0, 0, h, w);
                dI[1] = access(dIb, s1, 1, h, w);
                dI[2] = access(dIb, s2, 2, h, w);
                dI[3] = access(dIb, s3, 3, h, w);
                dI[4] = access(dIb, s4, 4, h, w);

                d0 = *(M+s0);
                d1 = 1/r*dx[1] + lambda*dI[1] + *(M+s1);
                d2 = 1/r*dx[2] + lambda*dI[2] + *(M+s2);
                d3 = 1/r*dx[3] + lambda*dI[3] + *(M+s3);
                d4 = 1/r*dx[4] + lambda*dI[4] + *(M+s4);

                double val[5] = {d0, d1, d2, d3, d4};
                map<int, double> idval;
                idval = minIDandVal(val, 5);
                map<int, double>::iterator it = idval.begin();
                switch((*it).first)
                {
                case 0: break;
                case 1: *(M+s0) = d1; *(L+s0) = *(L+s1); break;
                case 2: *(M+s0) = d2; *(L+s0) = *(L+s2); break;
                case 3: *(M+s0) = d3; *(L+s0) = *(L+s3); break;
                case 4: *(M+s0) = d4; *(L+s0) = *(L+s4); break;
                }

            }

            // case i = 1
            i = 0;
            s0 = i*h + j;
            s1 = (i+vbx[1])*h + j + vby[1];
            s2 = (i+vbx[2])*h + j + vby[2];
            s4 = (i+vbx[4])*h + j + vby[4];

            dI[0] = access(dIb, s0, 0, h, w);
            dI[1] = access(dIb, s1, 1, h, w);
            dI[2] = access(dIb, s2, 2, h, w);
            dI[4] = access(dIb, s3, 3, h, w);

            d0 = *(M+s0);
            d1 = 1/r*dx[1] + lambda*dI[1] + *(M+s1);
            d2 = 1/r*dx[2] + lambda*dI[2] + *(M+s2);
            d4 = 1/r*dx[4] + lambda*dI[4] + *(M+s4);

            val[0] = d0; val[1] =  d1; val[2] = d2; val[3] =  d4;
            idval = minIDandVal(val, 4);
            it = idval.begin();
            switch((*it).first)
            {
            case 0: break;
            case 1: *(M+s0) = d1; *(L+s0) = *(L+s1); break;
            case 2: *(M+s0) = d2; *(L+s0) = *(L+s2); break;
            case 3: *(M+s0) = d4; *(L+s0) = *(L+s4); break;
            }
        }



    }

    //contain values
    for(i=0;i<h*w;i++)
    {
    *(Mout+i) = *(M+i);
    *(Lout+i) = *(L+i);
    }

}
