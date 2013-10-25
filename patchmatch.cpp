/* ---------------------------------------------------------------------------
** This software is furnished "as is", without technical support,
** and with no warranty, express or implied, as to its usefulness for
** any purpose.
**
** patchmatch.cpp
** Matches patches between images
**
** Author: Sk. Mohammadul Haque
** Copyright (c) 2013 Sk. Mohammadul Haque
** For more details and updates, visit http://mohammadulhaque.alotspace.com
** -------------------------------------------------------------------------*/

#include <iostream>
#include <cmath>
#include <mex.h>
#include <vector>
#include <algorithm>
#include <stdint.h>
#include <climits>
#ifdef __PARALLEL__
#include <omp.h>
#endif
using namespace std;

#define DATA_TYPE int64_t
#define LNEG_VAL (INT_MIN)
#define LPOS_VAL (INT_MAX)

#define IS_REAL_2D_FULL_DOUBLE(P) (!mxIsComplex(P) && \
mxGetNumberOfDimensions(P)==2 && !mxIsSparse(P) && mxIsDouble(P))
#define IS_REAL_2D_OR_3D_FULL_DOUBLE(P) (!mxIsComplex(P) && \
(mxGetNumberOfDimensions(P)==2 || mxGetNumberOfDimensions(P)==3) && !mxIsEmpty(P) && !mxIsSparse(P) && mxIsDouble(P))
#define IS_REAL_SCALAR(P) (IS_REAL_2D_FULL_DOUBLE(P) && mxGetNumberOfElements(P)==1)
#define NNF_OUT plhs[0]
#define IMG_SRC_IN prhs[0]
#define IMG_DST_IN prhs[1]
#define PATCHSIZE_IN prhs[2]
#define NNK_IN prhs[3]
#define SEARCHRADIUS_IN prhs[4]
#define MASK_IN prhs[5]
#define SEARCHSTEP_IN prhs[6]
#define INCLUDESELF_IN prhs[7]
#define INCOMPLETE_IN prhs[8]
#define THRESHOLD_IN prhs[9]
#define DISTTYPE_IN prhs[10]

struct pxpair
{
    DATA_TYPE x, y;
};

#ifndef isnan
#define isnan(x) ((x)!=(x))
#endif

#ifndef min
#define min(x, y) ((x)<(y)?(x):(y))
#endif

typedef pxpair pxpair;

typedef vector<DATA_TYPE>::const_iterator vecintiter;
typedef vector<double>::const_iterator vecdoubleiter;
struct ordering
{
    bool operator ()(pair<DATA_TYPE, vecdoubleiter> const& a, pair<DATA_TYPE, vecdoubleiter> const& b)
    {
        return *(a.second) < *(b.second);
    }
};

template <typename T>vector<T> sort_from_ref( vector<T> const& in, vector<pair<DATA_TYPE, vecintiter> > const& reference)
{
    vector<T> ret(in.size());
    DATA_TYPE const size = in.size();
    for (DATA_TYPE i = 0; i < size; ++i)ret[i] = in[reference[i].first];
    return ret;
}

template <typename T>vector<T> sort_from_ref( vector<T> const& in, vector<pair<DATA_TYPE, vecdoubleiter> > const& reference)
{
    vector<T> ret(in.size());
    DATA_TYPE const size = in.size();
    for(DATA_TYPE i = 0; i < size; ++i) ret[i] = in[reference[i].first];
    return ret;
}

__inline double blockmatchgray(DATA_TYPE si, DATA_TYPE sj, DATA_TYPE ti, DATA_TYPE tj, double **mask, double *img_src, double *img_dst, DATA_TYPE patchsize, DATA_TYPE imgsM, DATA_TYPE imgsN, DATA_TYPE imgtM, DATA_TYPE imgtN);
__inline double blockmatchcolor(DATA_TYPE si, DATA_TYPE sj, DATA_TYPE ti, DATA_TYPE tj, double **mask, double *img_src, double *img_dst, DATA_TYPE patchsize, DATA_TYPE imgsM, DATA_TYPE imgsN, DATA_TYPE imgtM, DATA_TYPE imgtN);
__inline double blockmatchgray_one(DATA_TYPE si, DATA_TYPE sj, DATA_TYPE ti, DATA_TYPE tj, double **mask, double *img_src, double *img_dst, DATA_TYPE patchsize, DATA_TYPE imgsM, DATA_TYPE imgsN, DATA_TYPE imgtM, DATA_TYPE imgtN);
__inline double blockmatchcolor_one(DATA_TYPE si, DATA_TYPE sj, DATA_TYPE ti, DATA_TYPE tj, double **mask, double *img_src, double *img_dst, DATA_TYPE patchsize, DATA_TYPE imgsM, DATA_TYPE imgsN, DATA_TYPE imgtM, DATA_TYPE imgtN);
__inline double blockmatchgray_poisson(DATA_TYPE si, DATA_TYPE sj, DATA_TYPE ti, DATA_TYPE tj, double **mask, double *img_src, double *img_dst, DATA_TYPE patchsize, DATA_TYPE imgsM, DATA_TYPE imgsN, DATA_TYPE imgtM, DATA_TYPE imgtN);
__inline double blockmatchcolor_poisson(DATA_TYPE si, DATA_TYPE sj, DATA_TYPE ti, DATA_TYPE tj, double **mask, double *img_src, double *img_dst, DATA_TYPE patchsize, DATA_TYPE imgsM, DATA_TYPE imgsN, DATA_TYPE imgtM, DATA_TYPE imgtN);

__inline double blockmatchgray_cmu(DATA_TYPE si, DATA_TYPE sj, DATA_TYPE ti, DATA_TYPE tj, double **mask, double *img_src, double *img_dst, DATA_TYPE patchsize, DATA_TYPE imgsM, DATA_TYPE imgsN, DATA_TYPE imgtM, DATA_TYPE imgtN);
__inline double blockmatchcolor_cmu(DATA_TYPE si, DATA_TYPE sj, DATA_TYPE ti, DATA_TYPE tj, double **mask, double *img_src, double *img_dst, DATA_TYPE patchsize, DATA_TYPE imgsM, DATA_TYPE imgsN, DATA_TYPE imgtM, DATA_TYPE imgtN);

__inline double blockmatchgrayi(DATA_TYPE si, DATA_TYPE sj, DATA_TYPE ti, DATA_TYPE tj, double **mask, double *img, DATA_TYPE patchsize, DATA_TYPE imgM, DATA_TYPE imgN, double threshold);
__inline double blockmatchcolori(DATA_TYPE si, DATA_TYPE sj, DATA_TYPE ti, DATA_TYPE tj, double **mask, double *img, DATA_TYPE patchsize, DATA_TYPE imgM, DATA_TYPE imgN, double threshold);
__inline double blockmatchgrayi_one(DATA_TYPE si, DATA_TYPE sj, DATA_TYPE ti, DATA_TYPE tj, double **mask, double *img, DATA_TYPE patchsize, DATA_TYPE imgM, DATA_TYPE imgN, double threshold);
__inline double blockmatchcolori_one(DATA_TYPE si, DATA_TYPE sj, DATA_TYPE ti, DATA_TYPE tj, double **mask, double *img, DATA_TYPE patchsize, DATA_TYPE imgM, DATA_TYPE imgN, double threshold);
__inline double blockmatchgrayi_poisson(DATA_TYPE si, DATA_TYPE sj, DATA_TYPE ti, DATA_TYPE tj, double **mask, double *img, DATA_TYPE patchsize, DATA_TYPE imgM, DATA_TYPE imgN, double threshold);
__inline double blockmatchcolori_poisson(DATA_TYPE si, DATA_TYPE sj, DATA_TYPE ti, DATA_TYPE tj, double **mask, double *img, DATA_TYPE patchsize, DATA_TYPE imgM, DATA_TYPE imgN, double threshold);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double **mask = NULL, *rmask = NULL, *img_src = NULL, *img_dst = NULL, threshold = 0.0;
    int32_t *nnf = NULL;
    bool includeself = true, incomplete = false;
    int ndims_src, ndims_dst, nnk = 5, disttype = 2;
    DATA_TYPE IMG_DIMS_SRC[3]= {0, 0, 0}, IMG_DIMS_DST[3]= {0, 0, 0}, patchsize = 5, searchradius = 10, searchstep = 2;
    const mwSize *IMG_DIMS_SRC_ORI = NULL, *IMG_DIMS_DST_ORI = NULL;

    /* check number of arguments */
    /** nnf = patchmatch(img_src, img_dst, patchsize, nnk, searchradius, mask, searchstep, includeself, incomplete, threshold, disttype) **/

    if(nrhs<1 || nrhs>11) mexErrMsgTxt("Wrong number of input arguments.");
    else if(nlhs>1) mexErrMsgTxt("Too many output arguments.");

    /* check all arguments one by one if valid or set default if empty */

    if(!IS_REAL_2D_OR_3D_FULL_DOUBLE(IMG_SRC_IN)||mxIsEmpty(IMG_SRC_IN)) mexErrMsgTxt("IMG_SRC must be a real 2D or 3D full array.");
    if(!IS_REAL_2D_OR_3D_FULL_DOUBLE(IMG_DST_IN)&& !mxIsEmpty(IMG_DST_IN)) mexErrMsgTxt("IMG_DST must be a real 2D or 3D full array.");

    if(nrhs<11) disttype = 2;
    else
    {
        if(mxIsEmpty(DISTTYPE_IN)) disttype = 2;
        else
        {
            if(!IS_REAL_SCALAR(DISTTYPE_IN)) mexErrMsgTxt("DISTTYPE must be 0, 1, 2 or 3.");
            disttype = static_cast<int>(mxGetScalar(DISTTYPE_IN));
        }
    }

    if(nrhs<10) threshold = 0.0;
    else
    {
        if(mxIsEmpty(THRESHOLD_IN)) threshold = false;
        else
        {
            if(!IS_REAL_SCALAR(THRESHOLD_IN)) mexErrMsgTxt("THRESHOLD must be 0 and sum of mask.");
            threshold = static_cast<double>(mxGetScalar(THRESHOLD_IN));
        }
    }

    if(nrhs<9) incomplete = false;
    else
    {
        if(mxIsEmpty(INCOMPLETE_IN)) incomplete = false;
        else
        {
            if(!IS_REAL_SCALAR(INCOMPLETE_IN)) mexErrMsgTxt("INCOMPLETE must be 0 or 1.");
            incomplete = static_cast<bool>(mxGetScalar(INCOMPLETE_IN));
        }
    }

    if(nrhs<8) includeself = true;
    else
    {
        if(mxIsEmpty(INCLUDESELF_IN)) includeself = true;
        else
        {
            if(!IS_REAL_SCALAR(INCLUDESELF_IN)) mexErrMsgTxt("INCLUDESELF must be 0 or 1.");
            includeself = static_cast<bool>(mxGetScalar(INCLUDESELF_IN));
        }
    }

    if(nrhs<7) searchstep = 2;
    else
    {
        if(mxIsEmpty(SEARCHSTEP_IN)) searchstep = 2;
        else
        {
            if(!IS_REAL_SCALAR(SEARCHSTEP_IN)|| (mxGetScalar(SEARCHSTEP_IN)<1.0)) mexErrMsgTxt("SEARCHSTEP must be real positive scalar integer.");
            searchstep = static_cast<DATA_TYPE>(mxGetScalar(SEARCHSTEP_IN));
        }
    }

    if(nrhs<5) searchradius = 10;
    else
    {
        if(mxIsEmpty(SEARCHRADIUS_IN)) searchradius = 10;
        else
        {
            if(!IS_REAL_SCALAR(SEARCHRADIUS_IN)||mxGetScalar(SEARCHRADIUS_IN)<1.0) mexErrMsgTxt("SEARCHRADIUS must be real positive scalar integer or empty.");
            searchradius = static_cast<DATA_TYPE>(mxGetScalar(SEARCHRADIUS_IN));
        }
    }

    if(nrhs<4) nnk = 5;
    else
    {
        if(mxIsEmpty(NNK_IN)) nnk = 5;
        else
        {
            if(!IS_REAL_SCALAR(NNK_IN)||mxGetScalar(NNK_IN)<1.0) mexErrMsgTxt("NNK must be real positive scalar integer or empty.");
            nnk = static_cast<int>(mxGetScalar(NNK_IN));
        }
    }

    if(nrhs<3) patchsize = 5;
    else
    {
        if(mxIsEmpty(PATCHSIZE_IN)) patchsize = 5;
        else
        {
            if(!IS_REAL_SCALAR(PATCHSIZE_IN)||mxGetScalar(PATCHSIZE_IN)<1.0||(static_cast<int>(mxGetScalar(PATCHSIZE_IN))%2==0)) mexErrMsgTxt("PATCHSIZE must be real positive scalar odd integer or empty.");
            patchsize = mxGetScalar(PATCHSIZE_IN);
        }
    }

    /* get pointer to mask and its dimensions */
    /* ndims_src = mxGetNumberOfDimensions(MASK_IN); */ /* dont know why I added these unnecesssary lines */
    /* ndims_dst = mxGetNumberOfDimensions(MASK_IN); */ /* dont know why I added these unnecesssary lines */
    if(nrhs<6||mxIsEmpty(MASK_IN))
    {
        mask = new double*[patchsize];
        for(int i=0; i<patchsize; ++i)
        {
            mask[i] = new double[patchsize];
            for(int j=0; j<patchsize; ++j) mask[i][j] = 1.0;
        }
    }
    else
    {
        if(!IS_REAL_2D_FULL_DOUBLE(MASK_IN)) mexErrMsgTxt("MASK must be real full array or empty.");
        rmask = mxGetPr(MASK_IN);
        mask = new double*[patchsize];
        for(int i=0; i<patchsize; ++i)
        {
            mask[i] = new double[patchsize];
            for(int j=0; j<patchsize; ++j) mask[i][j] = rmask[i*patchsize+j];
        }
    }

    /* get pointer to src image and its dimensions */
    img_src = mxGetPr(IMG_SRC_IN);
    IMG_DIMS_SRC_ORI = mxGetDimensions(IMG_SRC_IN);
    ndims_src = mxGetNumberOfDimensions(IMG_SRC_IN);
    switch(ndims_src)
    {
    case 2:
        IMG_DIMS_SRC[0] = static_cast<DATA_TYPE>(IMG_DIMS_SRC_ORI[0]);
        IMG_DIMS_SRC[1] = static_cast<DATA_TYPE>(IMG_DIMS_SRC_ORI[1]);
        IMG_DIMS_SRC[2] = 1;
        break;
    case 3:
        IMG_DIMS_SRC[0] = static_cast<DATA_TYPE>(IMG_DIMS_SRC_ORI[0]);
        IMG_DIMS_SRC[1] = static_cast<DATA_TYPE>(IMG_DIMS_SRC_ORI[1]);
        IMG_DIMS_SRC[2] = static_cast<DATA_TYPE>(IMG_DIMS_SRC_ORI[2]);
        break;
    default:
        mexErrMsgTxt("IMG_SRC size error.");
        break;
    }

    if(patchsize>min(IMG_DIMS_SRC[0], IMG_DIMS_SRC[1])) mexErrMsgTxt("PATCHSIZE size is greater than IMG_SRC size.");

    /* get pointer to image and its dimensions */
    if(incomplete || mxIsEmpty(IMG_DST_IN))
    {
        img_dst = mxGetPr(IMG_SRC_IN);
        IMG_DIMS_DST_ORI = mxGetDimensions(IMG_SRC_IN);
        ndims_dst = mxGetNumberOfDimensions(IMG_SRC_IN);
    }
    else
    {
        img_dst = mxGetPr(IMG_DST_IN);
        IMG_DIMS_DST_ORI = mxGetDimensions(IMG_DST_IN);
        ndims_dst = mxGetNumberOfDimensions(IMG_DST_IN);
    }

    switch(ndims_dst)
    {
    case 2:
        IMG_DIMS_DST[0] = static_cast<DATA_TYPE>(IMG_DIMS_DST_ORI[0]);
        IMG_DIMS_DST[1] = static_cast<DATA_TYPE>(IMG_DIMS_DST_ORI[1]);
        IMG_DIMS_DST[2] = 1;
        break;
    case 3:
        IMG_DIMS_DST[0] = static_cast<DATA_TYPE>(IMG_DIMS_DST_ORI[0]);
        IMG_DIMS_DST[1] = static_cast<DATA_TYPE>(IMG_DIMS_DST_ORI[1]);
        IMG_DIMS_DST[2] = static_cast<DATA_TYPE>(IMG_DIMS_DST_ORI[2]);
        break;
    default:
        mexErrMsgTxt("IMG_DST size error.");
        break;
    }

    if(patchsize>min(IMG_DIMS_DST[0], IMG_DIMS_DST[1])) mexErrMsgTxt("PATCHSIZE size is greater than IMG_DST size.");

    mexPrintf("Src Img size : %" FMT64 "d %" FMT64 "d %" FMT64 "d\n", IMG_DIMS_SRC[0], IMG_DIMS_SRC[1], IMG_DIMS_SRC[2]);
    mexPrintf("Dst Img size : %" FMT64 "d %" FMT64 "d %" FMT64 "d\n", IMG_DIMS_DST[0], IMG_DIMS_DST[1], IMG_DIMS_DST[2]);
    mexPrintf("Patch size : %" FMT64 "d, Step size: %" FMT64 "d, Search radius: %" FMT64 "d\n", patchsize, searchstep, searchradius);
    if(IMG_DIMS_SRC[2]==1) mexPrintf("Type: Grayscale\n");
    if(IMG_DIMS_SRC[2]==3) mexPrintf("Type: Color\n");
    mexEvalString("drawnow;");

    DATA_TYPE max_X = IMG_DIMS_SRC[0]-patchsize;
    DATA_TYPE max_Y = IMG_DIMS_SRC[1]-patchsize;
    mwSize nnfdims[4] = {static_cast<mwSize>(IMG_DIMS_SRC[0]), static_cast<mwSize>(IMG_DIMS_SRC[1]), 3, static_cast<mwSize>(nnk)};
    NNF_OUT = mxCreateNumericArray(4, nnfdims, mxINT32_CLASS, mxREAL);
    nnf = reinterpret_cast<int32_t*>(mxGetPr(NNF_OUT));
    DATA_TYPE nelems = IMG_DIMS_SRC[0]*IMG_DIMS_SRC[1];

    if(!incomplete)
    {
        if(disttype == 3) /* poisson */
        {
            switch(IMG_DIMS_SRC[2])
            {
            case 1:
            {
#ifdef __PARALLEL__
                #pragma omp parallel for shared(IMG_DIMS_SRC, IMG_DIMS_DST, nnf, img_src, img_dst, searchradius, patchsize, searchstep, mask, nelems)
#endif

                for(DATA_TYPE sj=0; sj<=max_Y; ++sj)
                {
                    for(DATA_TYPE si=0; si<=max_X; ++si)
                    {
                        DATA_TYPE startx, starty, endx, endy;
                        DATA_TYPE Ms = IMG_DIMS_SRC[0], Ns = IMG_DIMS_SRC[1];
                        DATA_TYPE Mt = IMG_DIMS_DST[0], Nt = IMG_DIMS_DST[1];
                        vector<pxpair> pxall;
                        vector<double> distall;
                        double curr_dist;

                        startx = ((si*IMG_DIMS_DST[0])/IMG_DIMS_SRC[0])-searchradius;
                        startx = (startx>0)?startx: 0;
                        starty = ((sj*IMG_DIMS_DST[1])/IMG_DIMS_SRC[1])-searchradius;
                        starty = (starty>0)?starty: 0;

                        endx = ((si*IMG_DIMS_DST[0])/IMG_DIMS_SRC[0])+searchradius;
                        endx = (endx<(Mt-patchsize))?endx: (Mt-patchsize);
                        endy = ((sj*IMG_DIMS_DST[1])/IMG_DIMS_SRC[1])+searchradius;
                        endy = (endy<(Nt-patchsize))?endy: (Nt-patchsize);

                        for (DATA_TYPE nj = starty; nj<=endy; nj+=searchstep)
                        {
                            for (DATA_TYPE ni = startx; ni<=endx; ni+=searchstep)
                            {
                                if(si!=ni || sj!=nj||includeself)
                                {
                                    curr_dist = blockmatchgray_poisson(si, sj, ni, nj, mask, img_src, img_dst, patchsize, Ms, Ns, Mt, Nt);
                                    //mexPrintf("Source: (%Ld, %Ld) || Target: (%Ld, %Ld) || Dist: %f\n", si, sj, ni, nj, curr_dist);
                                    pxpair curr_px = {ni, nj};
                                    pxall.push_back(curr_px);
                                    distall.push_back(curr_dist);
                                }
                            }
                        }

                        //sort and take 1st few
                        vector<pair<DATA_TYPE, vecdoubleiter> > order(distall.size());

                        DATA_TYPE n = 0;
                        for (vecdoubleiter it = distall.begin(); it != distall.end(); ++it, ++n) order[n] = make_pair(n, it);
                        sort(order.begin(), order.end(), ordering());

                        vector<double> sorted_distall = sort_from_ref(distall, order);
                        vector<pxpair> sorted_pxall = sort_from_ref(pxall, order);
                        int lmax = (nnk<static_cast<int>(sorted_distall.size()))? nnk: static_cast<int>(sorted_distall.size());
                        for(int l=0; l<lmax; ++l)
                        {
                            nnf[si+sj*IMG_DIMS_SRC[0]+l*3*nelems] = sorted_pxall[l].y;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+1)*nelems] = sorted_pxall[l].x;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+2)*nelems] = static_cast<int32_t>(sorted_distall[l]);
                        }
                        for(int l=lmax; l<nnk; ++l)
                        {
                            nnf[si+sj*IMG_DIMS_SRC[0]+l*3*nelems] = LNEG_VAL;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+1)*nelems] = LNEG_VAL;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+2)*nelems] = LPOS_VAL;
                        }
                    }
                }
                break;
            }

            case 3:
            {
#ifdef __PARALLEL__
                #pragma omp parallel for shared(IMG_DIMS_SRC, IMG_DIMS_DST, nnf, img_src, img_dst, searchradius, patchsize, searchstep, mask, nelems)
#endif

                for(DATA_TYPE sj=0; sj<=max_Y; ++sj)
                {
                    for(DATA_TYPE si=0; si<=max_X; ++si)
                    {
                        DATA_TYPE startx, starty, endx, endy;
                        DATA_TYPE Ms = IMG_DIMS_SRC[0], Ns = IMG_DIMS_SRC[1];
                        DATA_TYPE Mt = IMG_DIMS_DST[0], Nt = IMG_DIMS_DST[1];
                        vector<pxpair> pxall;
                        vector<double> distall;
                        double curr_dist;

                        startx = ((si*IMG_DIMS_DST[0])/IMG_DIMS_SRC[0])-searchradius;
                        startx = (startx>0)?startx: 0;
                        starty = ((sj*IMG_DIMS_DST[1])/IMG_DIMS_SRC[1])-searchradius;
                        starty = (starty>0)?starty: 0;

                        endx = ((si*IMG_DIMS_DST[0])/IMG_DIMS_SRC[0])+searchradius;
                        endx = (endx<(Mt-patchsize))?endx: (Mt-patchsize);
                        endy = ((sj*IMG_DIMS_DST[1])/IMG_DIMS_SRC[1])+searchradius;
                        endy = (endy<(Nt-patchsize))?endy: (Nt-patchsize);

                        for (DATA_TYPE nj=starty; nj<=endy; nj+=searchstep)
                        {
                            for (DATA_TYPE ni=startx; ni<=endx; ni+=searchstep)
                            {
                                if(si!=ni || sj!=nj||includeself)
                                {
                                    curr_dist = blockmatchcolor_poisson(si, sj, ni, nj, mask, img_src, img_dst, patchsize, Ms, Ns, Mt, Nt);
                                    //mexPrintf("Source: (%Ld, %Ld) || Target: (%Ld, %Ld) || Dist: %f\n", si, sj, ni, nj, curr_dist);
                                    pxpair curr_px = {ni, nj};
                                    pxall.push_back(curr_px);
                                    distall.push_back(curr_dist);
                                }
                            }
                        }

                        //sort and take 1st few
                        vector<pair<DATA_TYPE, vecdoubleiter> > order(distall.size());

                        DATA_TYPE n = 0;
                        for (vecdoubleiter it = distall.begin(); it != distall.end(); ++it, ++n) order[n] = make_pair(n, it);
                        sort(order.begin(), order.end(), ordering());

                        vector<double> sorted_distall = sort_from_ref(distall, order);
                        vector<pxpair> sorted_pxall = sort_from_ref(pxall, order);
                        int lmax = (nnk<static_cast<int>(sorted_distall.size()))? nnk: static_cast<int>(sorted_distall.size());
                        for(int l=0; l<lmax; ++l)
                        {
                            nnf[si+sj*IMG_DIMS_SRC[0]+l*3*nelems] = sorted_pxall[l].y;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+1)*nelems] = sorted_pxall[l].x;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+2)*nelems] = static_cast<int32_t>(sorted_distall[l]);
                        }
                        for(int l=lmax; l<nnk; ++l)
                        {
                            nnf[si+sj*IMG_DIMS_SRC[0]+l*3*nelems] = LNEG_VAL;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+1)*nelems] = LNEG_VAL;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+2)*nelems] = LPOS_VAL;
                        }
                    }
                }
                break;
            }

            default:
            {
                // free mask memory created here before abnormal exit
                for(int i=0; i<patchsize; ++i)
                {
                    delete []mask[i];
                }
                delete []mask;
                mexErrMsgTxt("IMG 3rd dimension size must be 1 or 3.");
            }
            }
        }

        else if(disttype == 2) /* l2 */
        {
            switch(IMG_DIMS_SRC[2])
            {
            case 1:
            {
#ifdef __PARALLEL__
                #pragma omp parallel for shared(IMG_DIMS_SRC, IMG_DIMS_DST, nnf, img_src, img_dst, searchradius, patchsize, searchstep, mask, nelems)
#endif

                for(DATA_TYPE sj=0; sj<=max_Y; ++sj)
                {
                    for(DATA_TYPE si=0; si<=max_X; ++si)
                    {
                        DATA_TYPE startx, starty, endx, endy;
                        DATA_TYPE Ms = IMG_DIMS_SRC[0], Ns = IMG_DIMS_SRC[1];
                        DATA_TYPE Mt = IMG_DIMS_DST[0], Nt = IMG_DIMS_DST[1];
                        vector<pxpair> pxall;
                        vector<double> distall;
                        double curr_dist;

                        startx = ((si*IMG_DIMS_DST[0])/IMG_DIMS_SRC[0])-searchradius;
                        startx = (startx>0)?startx: 0;
                        starty = ((sj*IMG_DIMS_DST[1])/IMG_DIMS_SRC[1])-searchradius;
                        starty = (starty>0)?starty: 0;

                        endx = ((si*IMG_DIMS_DST[0])/IMG_DIMS_SRC[0])+searchradius;
                        endx = (endx<(Mt-patchsize))?endx: (Mt-patchsize);
                        endy = ((sj*IMG_DIMS_DST[1])/IMG_DIMS_SRC[1])+searchradius;
                        endy = (endy<(Nt-patchsize))?endy: (Nt-patchsize);

                        for (DATA_TYPE nj=starty; nj<=endy; nj+=searchstep)
                        {
                            for (DATA_TYPE ni=startx; ni<=endx; ni+=searchstep)
                            {
                                if(si!=ni || sj!=nj||includeself)
                                {
                                    curr_dist = blockmatchgray(si, sj, ni, nj, mask, img_src, img_dst, patchsize, Ms, Ns, Mt, Nt);
                                    //mexPrintf("Source: (%Ld, %Ld) || Target: (%Ld, %Ld) || Dist: %f\n", si, sj, ni, nj, curr_dist);
                                    pxpair curr_px = {ni, nj};
                                    pxall.push_back(curr_px);
                                    distall.push_back(curr_dist);
                                }
                            }
                        }

                        //sort and take 1st few
                        vector<pair<DATA_TYPE, vecdoubleiter> > order(distall.size());

                        DATA_TYPE n = 0;
                        for (vecdoubleiter it = distall.begin(); it != distall.end(); ++it, ++n) order[n] = make_pair(n, it);
                        sort(order.begin(), order.end(), ordering());

                        vector<double> sorted_distall = sort_from_ref(distall, order);
                        vector<pxpair> sorted_pxall = sort_from_ref(pxall, order);
                        int lmax = (nnk<static_cast<int>(sorted_distall.size()))? nnk: static_cast<int>(sorted_distall.size());
                        for(int l=0; l<lmax; ++l)
                        {
                            nnf[si+sj*IMG_DIMS_SRC[0]+l*3*nelems] = sorted_pxall[l].y;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+1)*nelems] = sorted_pxall[l].x;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+2)*nelems] = static_cast<int32_t>(sorted_distall[l]);
                        }
                        for(int l=lmax; l<nnk; ++l)
                        {
                            nnf[si+sj*IMG_DIMS_SRC[0]+l*3*nelems] = LNEG_VAL;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+1)*nelems] = LNEG_VAL;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+2)*nelems] = LPOS_VAL;
                        }
                    }
                }
                break;
            }

            case 3:
            {
#ifdef __PARALLEL__
                #pragma omp parallel for shared(IMG_DIMS_SRC, IMG_DIMS_DST, nnf, img_src, img_dst, searchradius, patchsize, searchstep, mask, nelems)
#endif

                for(DATA_TYPE sj=0; sj<=max_Y; ++sj)
                {
                    for(DATA_TYPE si=0; si<=max_X; ++si)
                    {
                        DATA_TYPE startx, starty, endx, endy;
                        DATA_TYPE Ms = IMG_DIMS_SRC[0], Ns = IMG_DIMS_SRC[1];
                        DATA_TYPE Mt = IMG_DIMS_DST[0], Nt = IMG_DIMS_DST[1];
                        vector<pxpair> pxall;
                        vector<double> distall;
                        double curr_dist;

                        startx = ((si*IMG_DIMS_DST[0])/IMG_DIMS_SRC[0])-searchradius;
                        startx = (startx>0)?startx: 0;
                        starty = ((sj*IMG_DIMS_DST[1])/IMG_DIMS_SRC[1])-searchradius;
                        starty = (starty>0)?starty: 0;

                        endx = ((si*IMG_DIMS_DST[0])/IMG_DIMS_SRC[0])+searchradius;
                        endx = (endx<(Mt-patchsize))?endx: (Mt-patchsize);
                        endy = ((sj*IMG_DIMS_DST[1])/IMG_DIMS_SRC[1])+searchradius;
                        endy = (endy<(Nt-patchsize))?endy: (Nt-patchsize);

                        for (DATA_TYPE nj=starty; nj<=endy; nj+=searchstep)
                        {
                            for (DATA_TYPE ni=startx; ni<=endx; ni+=searchstep)
                            {
                                if(si!=ni || sj!=nj||includeself)
                                {
                                    curr_dist = blockmatchcolor(si, sj, ni, nj, mask, img_src, img_dst, patchsize, Ms, Ns, Mt, Nt);
                                    //mexPrintf("Source: (%Ld, %Ld) || Target: (%Ld, %Ld) || Dist: %f\n", si, sj, ni, nj, curr_dist);
                                    pxpair curr_px = {ni, nj};
                                    pxall.push_back(curr_px);
                                    distall.push_back(curr_dist);
                                }
                            }
                        }

                        //sort and take 1st few
                        vector<pair<DATA_TYPE, vecdoubleiter> > order(distall.size());

                        DATA_TYPE n = 0;
                        for (vecdoubleiter it = distall.begin(); it != distall.end(); ++it, ++n) order[n] = make_pair(n, it);
                        sort(order.begin(), order.end(), ordering());

                        vector<double> sorted_distall = sort_from_ref(distall, order);
                        vector<pxpair> sorted_pxall = sort_from_ref(pxall, order);
                        int lmax = (nnk<static_cast<int>(sorted_distall.size()))? nnk: static_cast<int>(sorted_distall.size());
                        for(int l=0; l<lmax; ++l)
                        {
                            nnf[si+sj*IMG_DIMS_SRC[0]+l*3*nelems] = sorted_pxall[l].y;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+1)*nelems] = sorted_pxall[l].x;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+2)*nelems] = static_cast<int32_t>(sorted_distall[l]);
                        }
                        for(int l=lmax; l<nnk; ++l)
                        {
                            nnf[si+sj*IMG_DIMS_SRC[0]+l*3*nelems] = LNEG_VAL;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+1)*nelems] = LNEG_VAL;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+2)*nelems] = LPOS_VAL;
                        }
                    }
                }
                break;
            }

            default:
            {
                // free mask memory created here before abnormal exit
                for(int i=0; i<patchsize; ++i)
                {
                    delete []mask[i];
                }
                delete []mask;
                mexErrMsgTxt("IMG 3rd dimension size must be 1 or 3.");
            }
            }
        }

        else if(disttype == 1) /* l1 */
        {
            switch(IMG_DIMS_SRC[2])
            {
            case 1:
            {
#ifdef __PARALLEL__
                #pragma omp parallel for shared(IMG_DIMS_SRC, IMG_DIMS_DST, nnf, img_src, img_dst, searchradius, patchsize, searchstep, mask, nelems)
#endif

                for(DATA_TYPE sj=0; sj<=max_Y; ++sj)
                {
                    for(DATA_TYPE si=0; si<=max_X; ++si)
                    {
                        DATA_TYPE startx, starty, endx, endy;
                        DATA_TYPE Ms = IMG_DIMS_SRC[0], Ns = IMG_DIMS_SRC[1];
                        DATA_TYPE Mt = IMG_DIMS_DST[0], Nt = IMG_DIMS_DST[1];
                        vector<pxpair> pxall;
                        vector<double> distall;
                        double curr_dist;

                        startx = ((si*IMG_DIMS_DST[0])/IMG_DIMS_SRC[0])-searchradius;
                        startx = (startx>0)?startx: 0;
                        starty = ((sj*IMG_DIMS_DST[1])/IMG_DIMS_SRC[1])-searchradius;
                        starty = (starty>0)?starty: 0;

                        endx = ((si*IMG_DIMS_DST[0])/IMG_DIMS_SRC[0])+searchradius;
                        endx = (endx<(Mt-patchsize))?endx: (Mt-patchsize);
                        endy = ((sj*IMG_DIMS_DST[1])/IMG_DIMS_SRC[1])+searchradius;
                        endy = (endy<(Nt-patchsize))?endy: (Nt-patchsize);

                        for (DATA_TYPE nj = starty; nj<=endy; nj+=searchstep)
                        {
                            for (DATA_TYPE ni = startx; ni<=endx; ni+=searchstep)
                            {
                                if(si!=ni || sj!=nj||includeself)
                                {
                                    curr_dist = blockmatchgray_one(si, sj, ni, nj, mask, img_src, img_dst, patchsize, Ms, Ns, Mt, Nt);
                                    //mexPrintf("Source: (%Ld, %Ld) || Target: (%Ld, %Ld) || Dist: %f\n", si, sj, ni, nj, curr_dist);
                                    pxpair curr_px = {ni, nj};
                                    pxall.push_back(curr_px);
                                    distall.push_back(curr_dist);
                                }
                            }
                        }

                        //sort and take 1st few
                        vector<pair<DATA_TYPE, vecdoubleiter> > order(distall.size());

                        DATA_TYPE n = 0;
                        for (vecdoubleiter it = distall.begin(); it != distall.end(); ++it, ++n) order[n] = make_pair(n, it);
                        sort(order.begin(), order.end(), ordering());

                        vector<double> sorted_distall = sort_from_ref(distall, order);
                        vector<pxpair> sorted_pxall = sort_from_ref(pxall, order);
                        int lmax = (nnk<static_cast<int>(sorted_distall.size()))? nnk: static_cast<int>(sorted_distall.size());
                        for(int l=0; l<lmax; ++l)
                        {
                            nnf[si+sj*IMG_DIMS_SRC[0]+l*3*nelems] = sorted_pxall[l].y;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+1)*nelems] = sorted_pxall[l].x;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+2)*nelems] = static_cast<int32_t>(sorted_distall[l]);
                        }
                        for(int l=lmax; l<nnk; ++l)
                        {
                            nnf[si+sj*IMG_DIMS_SRC[0]+l*3*nelems] = LNEG_VAL;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+1)*nelems] = LNEG_VAL;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+2)*nelems] = LPOS_VAL;
                        }
                    }
                }
                break;
            }

            case 3:
            {
#ifdef __PARALLEL__
                #pragma omp parallel for shared(IMG_DIMS_SRC, IMG_DIMS_DST, nnf, img_src, img_dst, searchradius, patchsize, searchstep, mask, nelems)
#endif

                for(DATA_TYPE sj=0; sj<=max_Y; ++sj)
                {
                    for(DATA_TYPE si=0; si<=max_X; ++si)
                    {
                        DATA_TYPE startx, starty, endx, endy;
                        DATA_TYPE Ms = IMG_DIMS_SRC[0], Ns = IMG_DIMS_SRC[1];
                        DATA_TYPE Mt = IMG_DIMS_DST[0], Nt = IMG_DIMS_DST[1];
                        vector<pxpair> pxall;
                        vector<double> distall;
                        double curr_dist;

                        startx = ((si*IMG_DIMS_DST[0])/IMG_DIMS_SRC[0])-searchradius;
                        startx = (startx>0)?startx: 0;
                        starty = ((sj*IMG_DIMS_DST[1])/IMG_DIMS_SRC[1])-searchradius;
                        starty = (starty>0)?starty: 0;

                        endx = ((si*IMG_DIMS_DST[0])/IMG_DIMS_SRC[0])+searchradius;
                        endx = (endx<(Mt-patchsize))?endx: (Mt-patchsize);
                        endy = ((sj*IMG_DIMS_DST[1])/IMG_DIMS_SRC[1])+searchradius;
                        endy = (endy<(Nt-patchsize))?endy: (Nt-patchsize);

                        for (DATA_TYPE nj=starty; nj<=endy; nj+=searchstep)
                        {
                            for (DATA_TYPE ni=startx; ni<=endx; ni+=searchstep)
                            {
                                if(si!=ni || sj!=nj||includeself)
                                {
                                    curr_dist = blockmatchcolor_one(si, sj, ni, nj, mask, img_src, img_dst, patchsize, Ms, Ns, Mt, Nt);
                                    //mexPrintf("Source: (%Ld, %Ld) || Target: (%Ld, %Ld) || Dist: %f\n", si, sj, ni, nj, curr_dist);
                                    pxpair curr_px = {ni, nj};
                                    pxall.push_back(curr_px);
                                    distall.push_back(curr_dist);
                                }
                            }
                        }

                        //sort and take 1st few
                        vector<pair<DATA_TYPE, vecdoubleiter> > order(distall.size());

                        DATA_TYPE n = 0;
                        for (vecdoubleiter it = distall.begin(); it != distall.end(); ++it, ++n) order[n] = make_pair(n, it);
                        sort(order.begin(), order.end(), ordering());

                        vector<double> sorted_distall = sort_from_ref(distall, order);
                        vector<pxpair> sorted_pxall = sort_from_ref(pxall, order);
                        int lmax = (nnk<static_cast<int>(sorted_distall.size()))? nnk: static_cast<int>(sorted_distall.size());
                        for(int l=0; l<lmax; ++l)
                        {
                            nnf[si+sj*IMG_DIMS_SRC[0]+l*3*nelems] = sorted_pxall[l].y;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+1)*nelems] = sorted_pxall[l].x;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+2)*nelems] = static_cast<int32_t>(sorted_distall[l]);
                        }
                        for(int l=lmax; l<nnk; ++l)
                        {
                            nnf[si+sj*IMG_DIMS_SRC[0]+l*3*nelems] = LNEG_VAL;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+1)*nelems] = LNEG_VAL;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+2)*nelems] = LPOS_VAL;
                        }
                    }
                }
                break;
            }

            default:
            {
                // free mask memory created here before abnormal exit
                for(int i=0; i<patchsize; ++i)
                {
                    delete []mask[i];
                }
                delete []mask;
                mexErrMsgTxt("IMG 3rd dimension size must be 1 or 3.");
            }
            }
        }

        else if(disttype == 0) /* robust */
        {
            switch(IMG_DIMS_SRC[2])
            {
            case 1:
            {
#ifdef __PARALLEL__
                #pragma omp parallel for shared(IMG_DIMS_SRC, IMG_DIMS_DST, nnf, img_src, img_dst, searchradius, patchsize, searchstep, mask, nelems)
#endif

                for(DATA_TYPE sj=0; sj<=max_Y; ++sj)
                {
                    for(DATA_TYPE si=0; si<=max_X; ++si)
                    {
                        DATA_TYPE startx, starty, endx, endy;
                        DATA_TYPE Ms = IMG_DIMS_SRC[0], Ns = IMG_DIMS_SRC[1];
                        DATA_TYPE Mt = IMG_DIMS_DST[0], Nt = IMG_DIMS_DST[1];
                        vector<pxpair> pxall;
                        vector<double> distall;
                        double curr_dist;

                        startx = ((si*IMG_DIMS_DST[0])/IMG_DIMS_SRC[0])-searchradius;
                        startx = (startx>0)?startx: 0;
                        starty = ((sj*IMG_DIMS_DST[1])/IMG_DIMS_SRC[1])-searchradius;
                        starty = (starty>0)?starty: 0;

                        endx = ((si*IMG_DIMS_DST[0])/IMG_DIMS_SRC[0])+searchradius;
                        endx = (endx<(Mt-patchsize))?endx: (Mt-patchsize);
                        endy = ((sj*IMG_DIMS_DST[1])/IMG_DIMS_SRC[1])+searchradius;
                        endy = (endy<(Nt-patchsize))?endy: (Nt-patchsize);

                        for (DATA_TYPE nj = starty; nj<=endy; nj+=searchstep)
                        {
                            for (DATA_TYPE ni = startx; ni<=endx; ni+=searchstep)
                            {
                                if(si!=ni || sj!=nj||includeself)
                                {
                                    curr_dist = blockmatchgray_cmu(si, sj, ni, nj, mask, img_src, img_dst, patchsize, Ms, Ns, Mt, Nt);
                                    //mexPrintf("Source: (%Ld, %Ld) || Target: (%Ld, %Ld) || Dist: %f\n", si, sj, ni, nj, curr_dist);
                                    pxpair curr_px = {ni, nj};
                                    pxall.push_back(curr_px);
                                    distall.push_back(curr_dist);
                                }
                            }
                        }

                        //sort and take 1st few
                        vector<pair<DATA_TYPE, vecdoubleiter> > order(distall.size());

                        DATA_TYPE n = 0;
                        for (vecdoubleiter it = distall.begin(); it != distall.end(); ++it, ++n) order[n] = make_pair(n, it);
                        sort(order.begin(), order.end(), ordering());

                        vector<double> sorted_distall = sort_from_ref(distall, order);
                        vector<pxpair> sorted_pxall = sort_from_ref(pxall, order);
                        int lmax = (nnk<static_cast<int>(sorted_distall.size()))? nnk: static_cast<int>(sorted_distall.size());
                        for(int l=0; l<lmax; ++l)
                        {
                            nnf[si+sj*IMG_DIMS_SRC[0]+l*3*nelems] = sorted_pxall[l].y;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+1)*nelems] = sorted_pxall[l].x;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+2)*nelems] = static_cast<int32_t>(sorted_distall[l]);
                        }
                        for(int l=lmax; l<nnk; ++l)
                        {
                            nnf[si+sj*IMG_DIMS_SRC[0]+l*3*nelems] = LNEG_VAL;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+1)*nelems] = LNEG_VAL;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+2)*nelems] = LPOS_VAL;
                        }
                    }
                }
                break;
            }

            case 3:
            {
#ifdef __PARALLEL__
                #pragma omp parallel for shared(IMG_DIMS_SRC, IMG_DIMS_DST, nnf, img_src, img_dst, searchradius, patchsize, searchstep, mask, nelems)
#endif

                for(DATA_TYPE sj=0; sj<=max_Y; ++sj)
                {
                    for(DATA_TYPE si=0; si<=max_X; ++si)
                    {
                        DATA_TYPE startx, starty, endx, endy;
                        DATA_TYPE Ms = IMG_DIMS_SRC[0], Ns = IMG_DIMS_SRC[1];
                        DATA_TYPE Mt = IMG_DIMS_DST[0], Nt = IMG_DIMS_DST[1];
                        vector<pxpair> pxall;
                        vector<double> distall;
                        double curr_dist;

                        startx = ((si*IMG_DIMS_DST[0])/IMG_DIMS_SRC[0])-searchradius;
                        startx = (startx>0)?startx: 0;
                        starty = ((sj*IMG_DIMS_DST[1])/IMG_DIMS_SRC[1])-searchradius;
                        starty = (starty>0)?starty: 0;

                        endx = ((si*IMG_DIMS_DST[0])/IMG_DIMS_SRC[0])+searchradius;
                        endx = (endx<(Mt-patchsize))?endx: (Mt-patchsize);
                        endy = ((sj*IMG_DIMS_DST[1])/IMG_DIMS_SRC[1])+searchradius;
                        endy = (endy<(Nt-patchsize))?endy: (Nt-patchsize);

                        for (DATA_TYPE nj=starty; nj<=endy; nj+=searchstep)
                        {
                            for (DATA_TYPE ni=startx; ni<=endx; ni+=searchstep)
                            {
                                if(si!=ni || sj!=nj||includeself)
                                {
                                    curr_dist = blockmatchcolor_cmu(si, sj, ni, nj, mask, img_src, img_dst, patchsize, Ms, Ns, Mt, Nt);
                                    //mexPrintf("Source: (%Ld, %Ld) || Target: (%Ld, %Ld) || Dist: %f\n", si, sj, ni, nj, curr_dist);
                                    pxpair curr_px = {ni, nj};
                                    pxall.push_back(curr_px);
                                    distall.push_back(curr_dist);
                                }
                            }
                        }

                        //sort and take 1st few
                        vector<pair<DATA_TYPE, vecdoubleiter> > order(distall.size());

                        DATA_TYPE n = 0;
                        for (vecdoubleiter it = distall.begin(); it != distall.end(); ++it, ++n) order[n] = make_pair(n, it);
                        sort(order.begin(), order.end(), ordering());

                        vector<double> sorted_distall = sort_from_ref(distall, order);
                        vector<pxpair> sorted_pxall = sort_from_ref(pxall, order);
                        int lmax = (nnk<static_cast<int>(sorted_distall.size()))? nnk: static_cast<int>(sorted_distall.size());
                        for(int l=0; l<lmax; ++l)
                        {
                            nnf[si+sj*IMG_DIMS_SRC[0]+l*3*nelems] = sorted_pxall[l].y;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+1)*nelems] = sorted_pxall[l].x;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+2)*nelems] = static_cast<int32_t>(sorted_distall[l]);
                        }
                        for(int l=lmax; l<nnk; ++l)
                        {
                            nnf[si+sj*IMG_DIMS_SRC[0]+l*3*nelems] = LNEG_VAL;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+1)*nelems] = LNEG_VAL;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+2)*nelems] = LPOS_VAL;
                        }
                    }
                }
                break;
            }

            default:
            {
                // free mask memory created here before abnormal exit
                for(int i=0; i<patchsize; ++i)
                {
                    delete []mask[i];
                }
                delete []mask;
                mexErrMsgTxt("IMG 3rd dimension size must be 1 or 3.");
            }
            }
        }

        else
        {
            // free mask created here before exit
            for(int i=0; i<patchsize; ++i)
            {
                delete [] mask[i];
            }
            delete []mask;
            mexErrMsgTxt("Bad DISTTYPE given.");
        }


#ifdef __PARALLEL__
        #pragma omp parallel for shared(IMG_DIMS_SRC, nnf, nelems)
#endif
        for(DATA_TYPE sj=max_Y+1; sj<IMG_DIMS_SRC[1]; ++sj)
        {
            for(DATA_TYPE si=0; si<IMG_DIMS_SRC[0]; ++si)
            {
                for (DATA_TYPE l=0; l<nnk; ++l)
                {
                    nnf[si+sj*IMG_DIMS_SRC[0]+l*3*nelems] = LNEG_VAL;
                    nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+1)*nelems] = LNEG_VAL;
                    nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+2)*nelems] = LPOS_VAL;
                }
            }
        }

#ifdef __PARALLEL__
        #pragma omp parallel for shared(IMG_DIMS_SRC, nnf, nelems)
#endif
        for(DATA_TYPE sj=0; sj<=max_Y; ++sj)
        {
            for(DATA_TYPE si=max_X+1; si<IMG_DIMS_SRC[0]; ++si)
            {
                for (DATA_TYPE l=0; l<nnk; ++l)
                {
                    nnf[si+sj*IMG_DIMS_SRC[0]+l*3*nelems] = LNEG_VAL;
                    nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+1)*nelems] = LNEG_VAL;
                    nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+2)*nelems] = LPOS_VAL;
                }
            }
        }

    }
    /* nan / incomplete */ /* not done for two different images */
    else
    {
        if(disttype == 3) /* poisson */
        {
            switch(IMG_DIMS_SRC[2])
            {
            case 1:
            {
#ifdef __PARALLEL__
                #pragma omp parallel for shared(IMG_DIMS_SRC, IMG_DIMS_DST, nnf, img_src, img_dst, searchradius, patchsize, searchstep, mask, nelems, threshold)
#endif

                for(DATA_TYPE sj=0; sj<=max_Y; ++sj)
                {
                    for(DATA_TYPE si=0; si<=max_X; ++si)
                    {
                        DATA_TYPE startx, starty, endx, endy;
                        DATA_TYPE M = IMG_DIMS_SRC[0], N = IMG_DIMS_SRC[1];
                        vector<pxpair> pxall;
                        vector<double> distall;
                        double curr_dist, thresh = threshold;

                        startx = si-searchradius;
                        startx = (startx>0)?startx: 0;
                        starty = sj-searchradius;
                        starty = (starty>0)?starty: 0;

                        endx = si+searchradius;
                        endx = (endx<(M-patchsize))?endx: (M-patchsize);
                        endy = sj+searchradius;
                        endy = (endy<(N-patchsize))?endy: (N-patchsize);

                        for (DATA_TYPE nj=starty; nj<=endy; nj+=searchstep)
                        {
                            for (DATA_TYPE ni=startx; ni<=endx; ni+=searchstep)
                            {
                                if(si!=ni || sj!=nj||includeself)
                                {
                                    curr_dist = blockmatchgrayi_poisson(si, sj, ni, nj, mask, img_src, patchsize, M, N, thresh);
                                    //mexPrintf("Source: (%Ld, %Ld) || Target: (%Ld, %Ld) || Dist: %f\n", si, sj, ni, nj, curr_dist);
                                    pxpair curr_px = {ni, nj};
                                    pxall.push_back(curr_px);
                                    distall.push_back(curr_dist);
                                }
                            }
                        }

                        //sort and take 1st few
                        vector<pair<DATA_TYPE, vecdoubleiter> > order(distall.size());

                        DATA_TYPE n = 0;
                        for (vecdoubleiter it = distall.begin(); it != distall.end(); ++it, ++n) order[n] = make_pair(n, it);
                        sort(order.begin(), order.end(), ordering());

                        vector<double> sorted_distall = sort_from_ref(distall, order);
                        vector<pxpair> sorted_pxall = sort_from_ref(pxall, order);
                        int lmax = (nnk<static_cast<int>(sorted_distall.size()))? nnk: static_cast<int>(sorted_distall.size());
                        for(int l=0; l<lmax; ++l)
                        {
                            nnf[si+sj*IMG_DIMS_SRC[0]+l*3*nelems] = sorted_pxall[l].y;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+1)*nelems] = sorted_pxall[l].x;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+2)*nelems] = static_cast<int32_t>(sorted_distall[l]);
                        }
                        for(int l=lmax; l<nnk; ++l)
                        {
                            nnf[si+sj*IMG_DIMS_SRC[0]+l*3*nelems] = LNEG_VAL;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+1)*nelems] = LNEG_VAL;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+2)*nelems] = LPOS_VAL;
                        }
                    }
                }
                break;
            }

            case 3:
            {
#ifdef __PARALLEL__
                #pragma omp parallel for shared(IMG_DIMS_SRC, IMG_DIMS_DST, nnf, img_src, img_dst, searchradius, patchsize, searchstep, mask, nelems, threshold)
#endif

                for(DATA_TYPE sj=0; sj<=max_Y; ++sj)
                {
                    for(DATA_TYPE si=0; si<=max_X; ++si)
                    {
                        DATA_TYPE startx, starty, endx, endy;
                        DATA_TYPE M = IMG_DIMS_SRC[0], N = IMG_DIMS_SRC[1];
                        vector<pxpair> pxall;
                        vector<double> distall;
                        double curr_dist, thresh = threshold;

                        startx = si-searchradius;
                        startx = (startx>0)?startx: 0;
                        starty = sj-searchradius;
                        starty = (starty>0)?starty: 0;

                        endx = si+searchradius;
                        endx = (endx<(M-patchsize))?endx: (M-patchsize);
                        endy = sj+searchradius;
                        endy = (endy<(N-patchsize))?endy: (N-patchsize);

                        for (DATA_TYPE nj=starty; nj<=endy; nj+=searchstep)
                        {
                            for (DATA_TYPE ni=startx; ni<=endx; ni+=searchstep)
                            {
                                if(si!=ni || sj!=nj||includeself)
                                {
                                    curr_dist = blockmatchcolori_poisson(si, sj, ni, nj, mask, img_src, patchsize, M, N, thresh);
                                    //mexPrintf("Source: (%Ld, %Ld) || Target: (%Ld, %Ld) || Dist: %f\n", si, sj, ni, nj, curr_dist);
                                    pxpair curr_px = {ni, nj};
                                    pxall.push_back(curr_px);
                                    distall.push_back(curr_dist);
                                }
                            }
                        }

                        //sort and take 1st few
                        vector<pair<DATA_TYPE, vecdoubleiter> > order(distall.size());

                        DATA_TYPE n = 0;
                        for (vecdoubleiter it = distall.begin(); it != distall.end(); ++it, ++n) order[n] = make_pair(n, it);
                        sort(order.begin(), order.end(), ordering());

                        vector<double> sorted_distall = sort_from_ref(distall, order);
                        vector<pxpair> sorted_pxall = sort_from_ref(pxall, order);
                        int lmax = (nnk<static_cast<int>(sorted_distall.size()))? nnk: static_cast<int>(sorted_distall.size());
                        for(int l=0; l<lmax; ++l)
                        {
                            nnf[si+sj*IMG_DIMS_SRC[0]+l*3*nelems] = sorted_pxall[l].y;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+1)*nelems] = sorted_pxall[l].x;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+2)*nelems] = static_cast<int32_t>(sorted_distall[l]);
                        }
                        for(int l=lmax; l<nnk; ++l)
                        {
                            nnf[si+sj*IMG_DIMS_SRC[0]+l*3*nelems] = LNEG_VAL;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+1)*nelems] = LNEG_VAL;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+2)*nelems] = LPOS_VAL;
                        }
                    }
                }
                break;
            }

            default:
            {
                // free mask memory created here before abnormal exit
                for(int i=0; i<patchsize; ++i)
                {
                    delete []mask[i];
                }
                delete []mask;
                mexErrMsgTxt("IMG 3rd dimension size must be 1 or 3.");
            }
            }
        }

        else if(disttype == 2) /* l2 */
        {

            switch(IMG_DIMS_SRC[2])
            {
            case 1:
            {
#ifdef __PARALLEL__
                #pragma omp parallel for shared(IMG_DIMS_SRC, IMG_DIMS_DST, nnf, img_src, img_dst, searchradius, patchsize, searchstep, mask, nelems, threshold)
#endif

                for(DATA_TYPE sj=0; sj<=max_Y; ++sj)
                {
                    for(DATA_TYPE si=0; si<=max_X; ++si)
                    {
                        DATA_TYPE startx, starty, endx, endy;
                        DATA_TYPE M = IMG_DIMS_SRC[0], N = IMG_DIMS_SRC[1];
                        vector<pxpair> pxall;
                        vector<double> distall;
                        double curr_dist, thresh = threshold;

                        startx = si-searchradius;
                        startx = (startx>0)?startx: 0;
                        starty = sj-searchradius;
                        starty = (starty>0)?starty: 0;

                        endx = si+searchradius;
                        endx = (endx<(M-patchsize))?endx: (M-patchsize);
                        endy = sj+searchradius;
                        endy = (endy<(N-patchsize))?endy: (N-patchsize);

                        for (DATA_TYPE nj=starty; nj<=endy; nj+=searchstep)
                        {
                            for (DATA_TYPE ni=startx; ni<=endx; ni+=searchstep)
                            {
                                if(si!=ni || sj!=nj||includeself)
                                {
                                    curr_dist = blockmatchgrayi(si, sj, ni, nj, mask, img_src, patchsize, M, N, thresh);
                                    //mexPrintf("Source: (%Ld, %Ld) || Target: (%Ld, %Ld) || Dist: %f\n", si, sj, ni, nj, curr_dist);
                                    pxpair curr_px = {ni, nj};
                                    pxall.push_back(curr_px);
                                    distall.push_back(curr_dist);
                                }
                            }
                        }

                        //sort and take 1st few
                        vector<pair<DATA_TYPE, vecdoubleiter> > order(distall.size());

                        DATA_TYPE n = 0;
                        for (vecdoubleiter it = distall.begin(); it != distall.end(); ++it, ++n) order[n] = make_pair(n, it);
                        sort(order.begin(), order.end(), ordering());

                        vector<double> sorted_distall = sort_from_ref(distall, order);
                        vector<pxpair> sorted_pxall = sort_from_ref(pxall, order);
                        int lmax = (nnk<static_cast<int>(sorted_distall.size()))? nnk: static_cast<int>(sorted_distall.size());
                        for(int l=0; l<lmax; ++l)
                        {
                            nnf[si+sj*IMG_DIMS_SRC[0]+l*3*nelems] = sorted_pxall[l].y;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+1)*nelems] = sorted_pxall[l].x;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+2)*nelems] = static_cast<int32_t>(sorted_distall[l]);
                        }
                        for(int l=lmax; l<nnk; ++l)
                        {
                            nnf[si+sj*IMG_DIMS_SRC[0]+l*3*nelems] = LNEG_VAL;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+1)*nelems] = LNEG_VAL;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+2)*nelems] = LPOS_VAL;
                        }
                    }
                }
                break;
            }

            case 3:
            {
#ifdef __PARALLEL__
                #pragma omp parallel for shared(IMG_DIMS_SRC, IMG_DIMS_DST, nnf, img_src, img_dst, searchradius, patchsize, searchstep, mask, nelems, threshold)
#endif

                for(DATA_TYPE sj=0; sj<=max_Y; ++sj)
                {
                    for(DATA_TYPE si=0; si<=max_X; ++si)
                    {
                        DATA_TYPE startx, starty, endx, endy;
                        DATA_TYPE M = IMG_DIMS_SRC[0], N = IMG_DIMS_SRC[1];
                        vector<pxpair> pxall;
                        vector<double> distall;
                        double curr_dist, thresh = threshold;

                        startx = si-searchradius;
                        startx = (startx>0)?startx: 0;
                        starty = sj-searchradius;
                        starty = (starty>0)?starty: 0;

                        endx = si+searchradius;
                        endx = (endx<(M-patchsize))?endx: (M-patchsize);
                        endy = sj+searchradius;
                        endy = (endy<(N-patchsize))?endy: (N-patchsize);

                        for (DATA_TYPE nj=starty; nj<=endy; nj+=searchstep)
                        {
                            for (DATA_TYPE ni=startx; ni<=endx; ni+=searchstep)
                            {
                                if(si!=ni || sj!=nj||includeself)
                                {
                                    curr_dist = blockmatchcolori(si, sj, ni, nj, mask, img_src, patchsize, M, N, thresh);
                                    //mexPrintf("Source: (%Ld, %Ld) || Target: (%Ld, %Ld) || Dist: %f\n", si, sj, ni, nj, curr_dist);
                                    pxpair curr_px = {ni, nj};
                                    pxall.push_back(curr_px);
                                    distall.push_back(curr_dist);
                                }
                            }
                        }

                        //sort and take 1st few
                        vector<pair<DATA_TYPE, vecdoubleiter> > order(distall.size());

                        DATA_TYPE n = 0;
                        for (vecdoubleiter it = distall.begin(); it != distall.end(); ++it, ++n) order[n] = make_pair(n, it);
                        sort(order.begin(), order.end(), ordering());

                        vector<double> sorted_distall = sort_from_ref(distall, order);
                        vector<pxpair> sorted_pxall = sort_from_ref(pxall, order);
                        int lmax = (nnk<static_cast<int>(sorted_distall.size()))? nnk: static_cast<int>(sorted_distall.size());
                        for(int l=0; l<lmax; ++l)
                        {
                            nnf[si+sj*IMG_DIMS_SRC[0]+l*3*nelems] = sorted_pxall[l].y;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+1)*nelems] = sorted_pxall[l].x;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+2)*nelems] = static_cast<int32_t>(sorted_distall[l]);
                        }
                        for(int l=lmax; l<nnk; ++l)
                        {
                            nnf[si+sj*IMG_DIMS_SRC[0]+l*3*nelems] = LNEG_VAL;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+1)*nelems] = LNEG_VAL;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+2)*nelems] = LPOS_VAL;
                        }
                    }
                }
                break;
            }

            default:
            {
                // free mask memory created here before abnormal exit
                for(int i=0; i<patchsize; ++i)
                {
                    delete []mask[i];
                }
                delete []mask;
                mexErrMsgTxt("IMG 3rd dimension size must be 1 or 3.");
            }
            }
        }

        else if(disttype == 1) /* l1 */
        {
            switch(IMG_DIMS_SRC[2])
            {
            case 1:
            {
#ifdef __PARALLEL__
                #pragma omp parallel for shared(IMG_DIMS_SRC, IMG_DIMS_DST, nnf, img_src, img_dst, searchradius, patchsize, searchstep, mask, nelems, threshold)
#endif

                for(DATA_TYPE sj=0; sj<=max_Y; ++sj)
                {
                    for(DATA_TYPE si=0; si<=max_X; ++si)
                    {
                        DATA_TYPE startx, starty, endx, endy;
                        DATA_TYPE M = IMG_DIMS_SRC[0], N = IMG_DIMS_SRC[1];
                        vector<pxpair> pxall;
                        vector<double> distall;
                        double curr_dist, thresh = threshold;

                        startx = si-searchradius;
                        startx = (startx>0)?startx: 0;
                        starty = sj-searchradius;
                        starty = (starty>0)?starty: 0;

                        endx = si+searchradius;
                        endx = (endx<(M-patchsize))?endx: (M-patchsize);
                        endy = sj+searchradius;
                        endy = (endy<(N-patchsize))?endy: (N-patchsize);

                        for (DATA_TYPE nj=starty; nj<=endy; nj+=searchstep)
                        {
                            for (DATA_TYPE ni=startx; ni<=endx; ni+=searchstep)
                            {
                                if(si!=ni || sj!=nj||includeself)
                                {
                                    curr_dist = blockmatchgrayi_one(si, sj, ni, nj, mask, img_src, patchsize, M, N, thresh);
                                    //mexPrintf("Source: (%Ld, %Ld) || Target: (%Ld, %Ld) || Dist: %f\n", si, sj, ni, nj, curr_dist);
                                    pxpair curr_px = {ni, nj};
                                    pxall.push_back(curr_px);
                                    distall.push_back(curr_dist);
                                }
                            }
                        }

                        //sort and take 1st few
                        vector<pair<DATA_TYPE, vecdoubleiter> > order(distall.size());

                        DATA_TYPE n = 0;
                        for (vecdoubleiter it = distall.begin(); it != distall.end(); ++it, ++n) order[n] = make_pair(n, it);
                        sort(order.begin(), order.end(), ordering());

                        vector<double> sorted_distall = sort_from_ref(distall, order);
                        vector<pxpair> sorted_pxall = sort_from_ref(pxall, order);
                        int lmax = (nnk<static_cast<int>(sorted_distall.size()))? nnk: static_cast<int>(sorted_distall.size());
                        for(int l=0; l<lmax; ++l)
                        {
                            nnf[si+sj*IMG_DIMS_SRC[0]+l*3*nelems] = sorted_pxall[l].y;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+1)*nelems] = sorted_pxall[l].x;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+2)*nelems] = static_cast<int32_t>(sorted_distall[l]);
                        }
                        for(int l=lmax; l<nnk; ++l)
                        {
                            nnf[si+sj*IMG_DIMS_SRC[0]+l*3*nelems] = LNEG_VAL;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+1)*nelems] = LNEG_VAL;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+2)*nelems] = LPOS_VAL;
                        }
                    }
                }
                break;
            }

            case 3:
            {
#ifdef __PARALLEL__
                #pragma omp parallel for shared(IMG_DIMS_SRC, IMG_DIMS_DST, nnf, img_src, img_dst, searchradius, patchsize, searchstep, mask, nelems, threshold)
#endif

                for(DATA_TYPE sj=0; sj<=max_Y; ++sj)
                {
                    for(DATA_TYPE si=0; si<=max_X; ++si)
                    {
                        DATA_TYPE startx, starty, endx, endy;
                        DATA_TYPE M = IMG_DIMS_SRC[0], N = IMG_DIMS_SRC[1];
                        vector<pxpair> pxall;
                        vector<double> distall;
                        double curr_dist, thresh = threshold;

                        startx = si-searchradius;
                        startx = (startx>0)?startx: 0;
                        starty = sj-searchradius;
                        starty = (starty>0)?starty: 0;

                        endx = si+searchradius;
                        endx = (endx<(M-patchsize))?endx: (M-patchsize);
                        endy = sj+searchradius;
                        endy = (endy<(N-patchsize))?endy: (N-patchsize);

                        for (DATA_TYPE nj=starty; nj<=endy; nj+=searchstep)
                        {
                            for (DATA_TYPE ni=startx; ni<=endx; ni+=searchstep)
                            {
                                if(si!=ni || sj!=nj||includeself)
                                {
                                    curr_dist = blockmatchcolori_one(si, sj, ni, nj, mask, img_src, patchsize, M, N, thresh);
                                    //mexPrintf("Source: (%Ld, %Ld) || Target: (%Ld, %Ld) || Dist: %f\n", si, sj, ni, nj, curr_dist);
                                    pxpair curr_px = {ni, nj};
                                    pxall.push_back(curr_px);
                                    distall.push_back(curr_dist);
                                }
                            }
                        }

                        //sort and take 1st few
                        vector<pair<DATA_TYPE, vecdoubleiter> > order(distall.size());

                        DATA_TYPE n = 0;
                        for (vecdoubleiter it = distall.begin(); it != distall.end(); ++it, ++n) order[n] = make_pair(n, it);
                        sort(order.begin(), order.end(), ordering());

                        vector<double> sorted_distall = sort_from_ref(distall, order);
                        vector<pxpair> sorted_pxall = sort_from_ref(pxall, order);
                        int lmax = (nnk<static_cast<int>(sorted_distall.size()))? nnk: static_cast<int>(sorted_distall.size());
                        for(int l=0; l<lmax; ++l)
                        {
                            nnf[si+sj*IMG_DIMS_SRC[0]+l*3*nelems] = sorted_pxall[l].y;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+1)*nelems] = sorted_pxall[l].x;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+2)*nelems] = static_cast<int32_t>(sorted_distall[l]);
                        }
                        for(int l=lmax; l<nnk; ++l)
                        {
                            nnf[si+sj*IMG_DIMS_SRC[0]+l*3*nelems] = LNEG_VAL;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+1)*nelems] = LNEG_VAL;
                            nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+2)*nelems] = LPOS_VAL;
                        }
                    }
                }
                break;
            }

            default:
            {
                // free mask memory created here before abnormal exit
                for(int i=0; i<patchsize; ++i)
                {
                    delete []mask[i];
                }
                delete []mask;
                mexErrMsgTxt("IMG 3rd dimension size must be 1 or 3.");
            }
            }
        }

        else
        {
            // free mask created here before exit
            for(int i=0; i<patchsize; ++i)
            {
                delete [] mask[i];
            }
            delete []mask;
            mexErrMsgTxt("Bad DISTTYPE given.");
        }


#ifdef __PARALLEL__
        #pragma omp parallel for shared(IMG_DIMS_SRC, nnf, nelems)
#endif
        for(DATA_TYPE sj=max_Y+1; sj<IMG_DIMS_SRC[1]; ++sj)
        {
            for(DATA_TYPE si=0; si<IMG_DIMS_SRC[0]; ++si)
            {
                for (DATA_TYPE l=0; l<nnk; ++l)
                {
                    nnf[si+sj*IMG_DIMS_SRC[0]+l*3*nelems] = LNEG_VAL;
                    nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+1)*nelems] = LNEG_VAL;
                    nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+2)*nelems] = LPOS_VAL;
                }
            }
        }

#ifdef __PARALLEL__
        #pragma omp parallel for shared(IMG_DIMS_SRC, nnf, nelems)
#endif
        for(DATA_TYPE sj=0; sj<=max_Y; ++sj)
        {
            for(DATA_TYPE si=max_X+1; si<IMG_DIMS_SRC[0]; ++si)
            {
                for (DATA_TYPE l=0; l<nnk; ++l)
                {
                    nnf[si+sj*IMG_DIMS_SRC[0]+l*3*nelems] = LNEG_VAL;
                    nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+1)*nelems] = LNEG_VAL;
                    nnf[si+sj*IMG_DIMS_SRC[0]+(l*3+2)*nelems] = LPOS_VAL;
                }
            }
        }


    }

    // now free mask created here
    for(int i=0; i<patchsize; ++i)
    {
        delete [] mask[i];
    }
    delete []mask;
}


__inline double blockmatchgray(DATA_TYPE si, DATA_TYPE sj, DATA_TYPE ti, DATA_TYPE tj, double **mask, double *img_src, double *img_dst, DATA_TYPE patchsize, DATA_TYPE imgsM, DATA_TYPE imgsN, DATA_TYPE imgtM, DATA_TYPE imgtN)
{
    DATA_TYPE i, j, src_idx = sj*imgsM+si, tgt_idx = tj*imgtM+ti;
    double dist = 0.0;
    for(j=0; j<patchsize; ++j)
    {
        for(i=0; i<patchsize; ++i)
        {
            dist += mask[i][j]*(img_src[src_idx+i]-img_dst[tgt_idx+i])*(img_src[src_idx+i]-img_dst[tgt_idx+i]);
            //mexPrintf("(si, sj) = (%d, %d), (ti, tj) = (%d, %d),|| i = %d, j= %d, || sval = %f, tval = %f, dist = %f\n", si,sj,ti,tj, i,j,img[src_idx+i],img[tgt_idx+i], (img[src_idx+i]-img[tgt_idx+i])*(img[src_idx+i]-img[tgt_idx+i]));
        }
        src_idx += imgsM;
        tgt_idx += imgtM;
    }
    return dist;
}

__inline double blockmatchcolor(DATA_TYPE si, DATA_TYPE sj, DATA_TYPE ti, DATA_TYPE tj, double **mask, double *img_src, double *img_dst, DATA_TYPE patchsize, DATA_TYPE imgsM, DATA_TYPE imgsN, DATA_TYPE imgtM, DATA_TYPE imgtN)
{
    DATA_TYPE i, j, src_idx = sj*imgsM+si, tgt_idx = tj*imgtM+ti, numels_src = imgsM*imgsN, numels_dst = imgtM*imgtN;
    double dist = 0.0;
    /* channel R */
    for(j=0; j<patchsize; ++j)
    {
        for(i=0; i<patchsize; ++i)
        {
            dist += mask[i][j]*(img_src[src_idx+i]-img_dst[tgt_idx+i])*(img_src[src_idx+i]-img_dst[tgt_idx+i]);
        }
        src_idx += imgsM;
        tgt_idx += imgtM;
    }

    /* channel G */
    src_idx = sj*imgsM+si+numels_src;
    tgt_idx = tj*imgtM+ti+numels_dst; /* shift to the next layer */
    for(j=0; j<patchsize; ++j)
    {
        for(i=0; i<patchsize; ++i)
        {
            dist += mask[i][j]*(img_src[src_idx+i]-img_dst[tgt_idx+i])*(img_src[src_idx+i]-img_dst[tgt_idx+i]);
        }
        src_idx += imgsM;
        tgt_idx += imgtM;
    }

    /* channel B */
    src_idx = sj*imgsM+si+2*numels_src;
    tgt_idx = tj*imgtM+ti+2*numels_dst; /* shift to the next layer */
    for(j=0; j<patchsize; ++j)
    {
        for(i=0; i<patchsize; ++i)
        {
            dist += mask[i][j]*(img_src[src_idx+i]-img_dst[tgt_idx+i])*(img_src[src_idx+i]-img_dst[tgt_idx+i]);
        }
        src_idx += imgsM;
        tgt_idx += imgtM;
    }
    return dist;
}



__inline double blockmatchgray_one(DATA_TYPE si, DATA_TYPE sj, DATA_TYPE ti, DATA_TYPE tj, double **mask, double *img_src, double *img_dst, DATA_TYPE patchsize, DATA_TYPE imgsM, DATA_TYPE imgsN, DATA_TYPE imgtM, DATA_TYPE imgtN)
{
    DATA_TYPE i, j, src_idx = sj*imgsM+si, tgt_idx = tj*imgtM+ti;
    double dist = 0.0;
    for(j=0; j<patchsize; ++j)
    {
        for(i=0; i<patchsize; ++i)
        {
            dist += mask[i][j]*fabs(img_src[src_idx+i]-img_dst[tgt_idx+i]);
            //mexPrintf("(si, sj) = (%d, %d), (ti, tj) = (%d, %d),|| i = %d, j= %d, || sval = %f, tval = %f, dist = %f\n", si,sj,ti,tj, i,j,img[src_idx+i],img[tgt_idx+i], (img[src_idx+i]-img[tgt_idx+i])*(img[src_idx+i]-img[tgt_idx+i]));
        }
        src_idx += imgsM;
        tgt_idx += imgtM;
    }
    return dist;
}

__inline double blockmatchcolor_one(DATA_TYPE si, DATA_TYPE sj, DATA_TYPE ti, DATA_TYPE tj, double **mask, double *img_src, double *img_dst, DATA_TYPE patchsize, DATA_TYPE imgsM, DATA_TYPE imgsN, DATA_TYPE imgtM, DATA_TYPE imgtN)
{
    DATA_TYPE i, j, src_idx = sj*imgsM+si, tgt_idx = tj*imgtM+ti, numels_src = imgsM*imgsN, numels_dst = imgtM*imgtN;
    double dist = 0.0;
    /* channel R */
    for(j=0; j<patchsize; ++j)
    {
        for(i=0; i<patchsize; ++i)
        {
            dist += mask[i][j]*fabs(img_src[src_idx+i]-img_dst[tgt_idx+i]);
        }
        src_idx += imgsM;
        tgt_idx += imgtM;
    }

    /* channel G */
    src_idx = sj*imgsM+si+numels_src;
    tgt_idx = tj*imgtM+ti+numels_dst; /* shift to the next layer */
    for(j=0; j<patchsize; ++j)
    {
        for(i=0; i<patchsize; ++i)
        {
            dist += mask[i][j]*fabs(img_src[src_idx+i]-img_dst[tgt_idx+i]);
        }
        src_idx += imgsM;
        tgt_idx += imgtM;
    }

    /* channel B */
    src_idx = sj*imgsM+si+2*numels_src;
    tgt_idx = tj*imgtM+ti+2*numels_dst; /* shift to the next layer */
    for(j=0; j<patchsize; ++j)
    {
        for(i=0; i<patchsize; ++i)
        {
            dist += mask[i][j]*fabs(img_src[src_idx+i]-img_dst[tgt_idx+i]);
        }
        src_idx += imgsM;
        tgt_idx += imgtM;
    }
    return dist;
}


__inline double blockmatchgray_cmu(DATA_TYPE si, DATA_TYPE sj, DATA_TYPE ti, DATA_TYPE tj, double **mask, double *img_src, double *img_dst, DATA_TYPE patchsize, DATA_TYPE imgsM, DATA_TYPE imgsN, DATA_TYPE imgtM, DATA_TYPE imgtN)
{
    DATA_TYPE i, j, src_idx = sj*imgsM+si, p_m_one = patchsize-1;
    double dist = 0.0, tempxx, tempxy, tempyy;
    double s_C[3] = {0,0,0};

    //calculate Cs
    src_idx = sj*imgsM+si;

    for(j=0; j<p_m_one; ++j)
    {
        for(i=0; i<p_m_one; ++i)
        {
            s_C[0] += (img_src[src_idx+i+1]-img_src[src_idx+i])*(img_src[src_idx+i+1]-img_src[src_idx+i]); //along x
            s_C[1] += (img_src[src_idx+i+imgsM]-img_src[src_idx+i])*(img_src[src_idx+i+imgsM]-img_src[src_idx+i]); //along y
            s_C[2] += (img_src[src_idx+i+1]-img_src[src_idx+i])*(img_src[src_idx+i+imgsM]-img_src[src_idx+i]); //along x-y

        }

        src_idx += imgsM+1;
    }

    s_C[0] = (s_C[0])/(p_m_one*p_m_one);
    s_C[1] = (s_C[1])/(p_m_one*p_m_one);
    s_C[2] = (s_C[2])/(p_m_one*p_m_one);
    tempxx = static_cast<double>(si-ti); tempxx *= tempxx;
    tempyy = static_cast<double>(sj-tj); tempyy *= tempyy;
    tempxy = static_cast<double>(si-ti); tempxy *= static_cast<double>(sj-tj);


    dist = s_C[0]*tempxx+s_C[1]*tempyy+2*s_C[2]*tempxy;
    //mexPrintf("(si, sj) = (%d, %d), (ti, tj) = (%d, %d),|| i = %d, j= %d, || sval = %f, tval = %f, dist = %f\n", si,sj,ti,tj, i,j,img[src_idx+i],img[tgt_idx+i], (img[src_idx+i]-img[tgt_idx+i])*(img[src_idx+i]-img[tgt_idx+i]));
    return dist;
}

__inline double blockmatchcolor_cmu(DATA_TYPE si, DATA_TYPE sj, DATA_TYPE ti, DATA_TYPE tj, double **mask, double *img_src, double *img_dst, DATA_TYPE patchsize, DATA_TYPE imgsM, DATA_TYPE imgsN, DATA_TYPE imgtM, DATA_TYPE imgtN)
{
    DATA_TYPE i, j, src_idx = sj*imgsM+si, p_m_one = patchsize-1, numels = imgsM*imgsN;
    double dist = 0.0, tempxx, tempxy, tempyy;
    double s_C[3][3] = {{0,0,0},{0,0,0},{0,0,0}};

    //for R channel

    //calculate Cs
    for(j=0; j<p_m_one; ++j)
    {
        for(i=0; i<p_m_one; ++i)
        {
            s_C[0][0] += (img_src[src_idx+i+1]-img_src[src_idx+i])*(img_src[src_idx+i+1]-img_src[src_idx+i]); //along x
            s_C[0][1] += (img_src[src_idx+i+imgsM]-img_src[src_idx+i])*(img_src[src_idx+i+imgsM]-img_src[src_idx+i]); //along y
            s_C[0][2] += (img_src[src_idx+i+1]-img_src[src_idx+i])*(img_src[src_idx+i+imgsM]-img_src[src_idx+i]); //along x-y

        }
        src_idx += imgsM+1;
    }

    s_C[0][0] = (s_C[0][0])/(p_m_one*p_m_one);
    s_C[0][1] = (s_C[0][1])/(p_m_one*p_m_one);
    s_C[0][2] = (s_C[0][2])/(p_m_one*p_m_one);

    //for G channel
    src_idx = sj*imgsM+si+numels; /* shift to the next layer */

    //calculate Cs
    for(j=0; j<p_m_one; ++j)
    {
        for(i=0; i<p_m_one; ++i)
        {
            s_C[1][0] += (img_src[src_idx+i+1]-img_src[src_idx+i])*(img_src[src_idx+i+1]-img_src[src_idx+i]); //along x
            s_C[1][1] += (img_src[src_idx+i+imgsM]-img_src[src_idx+i])*(img_src[src_idx+i+imgsM]-img_src[src_idx+i]); //along y
            s_C[1][2] += (img_src[src_idx+i+1]-img_src[src_idx+i])*(img_src[src_idx+i+imgsM]-img_src[src_idx+i]); //along x-y

        }
        src_idx += imgsM+1;
    }

    s_C[1][0] = (s_C[1][0])/(p_m_one*p_m_one);
    s_C[1][1] = (s_C[1][1])/(p_m_one*p_m_one);
    s_C[1][2] = (s_C[1][2])/(p_m_one*p_m_one);

    //for B channel
    src_idx = sj*imgsM+si+2*numels; /* shift to the next layer */

    //calculate Cs
    for(j=0; j<p_m_one; ++j)
    {
        for(i=0; i<p_m_one; ++i)
        {
            s_C[2][0] += (img_src[src_idx+i+1]-img_src[src_idx+i])*(img_src[src_idx+i+1]-img_src[src_idx+i]); //along x
            s_C[2][1] += (img_src[src_idx+i+imgsM]-img_src[src_idx+i])*(img_src[src_idx+i+imgsM]-img_src[src_idx+i]); //along y
            s_C[2][2] += (img_src[src_idx+i+1]-img_src[src_idx+i])*(img_src[src_idx+i+imgsM]-img_src[src_idx+i]); //along x-y
        }
        src_idx += imgsM+1;
    }

    s_C[2][0] = (s_C[2][0])/(p_m_one*p_m_one);
    s_C[2][1] = (s_C[2][1])/(p_m_one*p_m_one);
    s_C[2][2] = (s_C[2][2])/(p_m_one*p_m_one);

    tempxx = static_cast<double>(si-ti); tempxx *= tempxx;
    tempyy = static_cast<double>(sj-tj); tempyy *= tempyy;
    tempxy = static_cast<double>(si-ti); tempxy *= static_cast<double>(sj-tj);

    dist = (s_C[0][0]+s_C[1][0]+s_C[2][0])*tempxx+(s_C[0][1]+s_C[1][1]+s_C[2][1])*tempyy+2*(s_C[0][2]+s_C[1][2]+s_C[2][2])*tempxy;
    //mexPrintf("(si, sj) = (%d, %d), (ti, tj) = (%d, %d),|| i = %d, j= %d, || sval = %f, tval = %f, dist = %f\n", si,sj,ti,tj, i,j,img[src_idx+i],img[tgt_idx+i], (img[src_idx+i]-img[tgt_idx+i])*(img[src_idx+i]-img[tgt_idx+i]));

    return dist;
}



__inline double blockmatchgray_poisson(DATA_TYPE si, DATA_TYPE sj, DATA_TYPE ti, DATA_TYPE tj, double **mask, double *img_src, double *img_dst, DATA_TYPE patchsize, DATA_TYPE imgsM, DATA_TYPE imgsN, DATA_TYPE imgtM, DATA_TYPE imgtN)
{
    DATA_TYPE i, j, src_idx = sj*imgsM+si, tgt_idx = tj*imgtM+ti;
    double dist = 0.0;
    for(j=0; j<patchsize; ++j)
    {
        for(i=0; i<patchsize; ++i)
        {
            dist += mask[i][j]*(lgamma(img_src[src_idx+i]+img_dst[tgt_idx+i]+0.5)/sqrt(lgamma(2*img_src[src_idx+i]+0.5)*lgamma(2*img_dst[tgt_idx+i]+0.5)));
            //mexPrintf("(si, sj) = (%d, %d), (ti, tj) = (%d, %d),|| i = %d, j= %d, || sval = %f, tval = %f, dist = %f\n", si,sj,ti,tj, i,j,img[src_idx+i],img[tgt_idx+i], (img[src_idx+i]-img[tgt_idx+i])*(img[src_idx+i]-img[tgt_idx+i]));
        }
        src_idx += imgsM;
        tgt_idx += imgtM;
    }
    return dist;
}

__inline double blockmatchcolor_poisson(DATA_TYPE si, DATA_TYPE sj, DATA_TYPE ti, DATA_TYPE tj, double **mask, double *img_src, double *img_dst, DATA_TYPE patchsize, DATA_TYPE imgsM, DATA_TYPE imgsN, DATA_TYPE imgtM, DATA_TYPE imgtN)
{
    DATA_TYPE i, j, src_idx = sj*imgsM+si, tgt_idx = tj*imgtM+ti, numels_src = imgsM*imgsN, numels_dst = imgtM*imgtN;
    double dist = 0.0;
    /* channel R */
    for(j=0; j<patchsize; ++j)
    {
        for(i=0; i<patchsize; ++i)
        {
            dist += mask[i][j]*(lgamma(img_src[src_idx+i]+img_dst[tgt_idx+i]+0.5)/sqrt(lgamma(2*img_src[src_idx+i]+0.5)*lgamma(2*img_dst[tgt_idx+i]+0.5)));
        }
        src_idx += imgsM;
        tgt_idx += imgtM;
    }

    /* channel G */
    src_idx = sj*imgsM+si+numels_src;
    tgt_idx = tj*imgtM+ti+numels_dst; /* shift to the next layer */
    for(j=0; j<patchsize; ++j)
    {
        for(i=0; i<patchsize; ++i)
        {
            dist += mask[i][j]*(lgamma(img_src[src_idx+i]+img_dst[tgt_idx+i]+0.5)/sqrt(lgamma(2*img_src[src_idx+i]+0.5)*lgamma(2*img_dst[tgt_idx+i]+0.5)));
        }
        src_idx += imgsM;
        tgt_idx += imgtM;
    }

    /* channel B */
    src_idx = sj*imgsM+si+2*numels_src;
    tgt_idx = tj*imgtM+ti+2*numels_dst; /* shift to the next layer */
    for(j=0; j<patchsize; ++j)
    {
        for(i=0; i<patchsize; ++i)
        {
            dist += mask[i][j]*(lgamma(img_src[src_idx+i]+img_dst[tgt_idx+i]+0.5)/sqrt(lgamma(2*img_src[src_idx+i]+0.5)*lgamma(2*img_dst[tgt_idx+i]+0.5)));
        }
        src_idx += imgsM;
        tgt_idx += imgtM;
    }
    return dist;
}




__inline double blockmatchgrayi(DATA_TYPE si, DATA_TYPE sj, DATA_TYPE ti, DATA_TYPE tj, double **mask, double *img, DATA_TYPE patchsize, DATA_TYPE imgM, DATA_TYPE imgN, double threshold)
{
    DATA_TYPE i, j, src_idx = sj*imgM+si, tgt_idx = tj*imgM+ti;
    double dist = 0.0, den = 0, temp;
    for(j=0; j<patchsize; ++j)
    {
        for(i=0; i<patchsize; ++i)
        {
            temp = (img[src_idx+i]-img[tgt_idx+i])*(img[src_idx+i]-img[tgt_idx+i]);
            if(!isnan(temp))
            {
                dist += mask[i][j]*temp;
                den += mask[i][j];
            }
            //mexPrintf("(si, sj) = (%d, %d), (ti, tj) = (%d, %d),|| i = %d, j= %d, || sval = %f, tval = %f, dist = %f\n", si,sj,ti,tj, i,j,img[src_idx+i],img[tgt_idx+i], (img[src_idx+i]-img[tgt_idx+i])*(img[src_idx+i]-img[tgt_idx+i]));
        }
        src_idx += imgM;
        tgt_idx += imgM;
    }
    dist /= den;
    if(den<=threshold) dist = static_cast<double>(LPOS_VAL);
    return dist;
}

__inline double blockmatchcolori(DATA_TYPE si, DATA_TYPE sj, DATA_TYPE ti, DATA_TYPE tj, double **mask, double *img, DATA_TYPE patchsize, DATA_TYPE imgM, DATA_TYPE imgN, double threshold)
{
    DATA_TYPE i, j, src_idx = sj*imgM+si, tgt_idx = tj*imgM+ti, numels = imgM*imgN;
    double dist = 0.0, den = 0, temp;
    /* channel R */
    for(j=0; j<patchsize; ++j)
    {
        for(i=0; i<patchsize; ++i)
        {
            temp = mask[i][j]*(img[src_idx+i]-img[tgt_idx+i])*(img[src_idx+i]-img[tgt_idx+i]);
            if(!isnan(temp))
            {
                dist += mask[i][j]*temp;
                den += mask[i][j];
            }
        }
        src_idx += imgM;
        tgt_idx += imgM;
    }

    /* channel G */
    src_idx = sj*imgM+si+numels;
    tgt_idx = tj*imgM+ti+numels; /* shift to the next layer */
    for(j=0; j<patchsize; ++j)
    {
        for(i=0; i<patchsize; ++i)
        {
            temp = mask[i][j]*(img[src_idx+i]-img[tgt_idx+i])*(img[src_idx+i]-img[tgt_idx+i]);
            if(!isnan(temp))
            {
                dist += mask[i][j]*temp;
                den += mask[i][j];
            }
        }
        src_idx += imgM;
        tgt_idx += imgM;
    }

    /* channel B */
    src_idx = sj*imgM+si+2*numels;
    tgt_idx = tj*imgM+ti+2*numels; /* shift to the next layer */
    for(j=0; j<patchsize; ++j)
    {
        for(i=0; i<patchsize; ++i)
        {
            temp = mask[i][j]*(img[src_idx+i]-img[tgt_idx+i])*(img[src_idx+i]-img[tgt_idx+i]);
            if(!isnan(temp))
            {
                dist += mask[i][j]*temp;
                den += mask[i][j];
            }
        }
        src_idx += imgM;
        tgt_idx += imgM;
    }
    dist /= den;
    if(den<=threshold) dist = static_cast<double>(LPOS_VAL);
    return dist;
}


__inline double blockmatchgrayi_one(DATA_TYPE si, DATA_TYPE sj, DATA_TYPE ti, DATA_TYPE tj, double **mask, double *img, DATA_TYPE patchsize, DATA_TYPE imgM, DATA_TYPE imgN, double threshold)
{
    DATA_TYPE i, j, src_idx = sj*imgM+si, tgt_idx = tj*imgM+ti;
    double dist = 0.0, den = 0, temp;
    for(j=0; j<patchsize; ++j)
    {
        for(i=0; i<patchsize; ++i)
        {
            temp = fabs(img[src_idx+i]-img[tgt_idx+i]);
            if(!isnan(temp))
            {
                dist += mask[i][j]*temp;
                den += mask[i][j];
            }
            //mexPrintf("(si, sj) = (%d, %d), (ti, tj) = (%d, %d),|| i = %d, j= %d, || sval = %f, tval = %f, dist = %f\n", si,sj,ti,tj, i,j,img[src_idx+i],img[tgt_idx+i], (img[src_idx+i]-img[tgt_idx+i])*(img[src_idx+i]-img[tgt_idx+i]));
        }
        src_idx += imgM;
        tgt_idx += imgM;
    }
    dist /= den;
    if(den<=threshold) dist = static_cast<double>(LPOS_VAL);
    return dist;
}

__inline double blockmatchcolori_one(DATA_TYPE si, DATA_TYPE sj, DATA_TYPE ti, DATA_TYPE tj, double **mask, double *img, DATA_TYPE patchsize, DATA_TYPE imgM, DATA_TYPE imgN, double threshold)
{
    DATA_TYPE i, j, src_idx = sj*imgM+si, tgt_idx = tj*imgM+ti, numels = imgM*imgN;
    double dist = 0.0, den = 0, temp;
    /* channel R */
    for(j=0; j<patchsize; ++j)
    {
        for(i=0; i<patchsize; ++i)
        {
            temp = mask[i][j]*fabs(img[src_idx+i]-img[tgt_idx+i]);
            if(!isnan(temp))
            {
                dist += mask[i][j]*temp;
                den += mask[i][j];
            }
        }
        src_idx += imgM;
        tgt_idx += imgM;
    }

    /* channel G */
    src_idx = sj*imgM+si+numels;
    tgt_idx = tj*imgM+ti+numels; /* shift to the next layer */
    for(j=0; j<patchsize; ++j)
    {
        for(i=0; i<patchsize; ++i)
        {
            temp = mask[i][j]*fabs(img[src_idx+i]-img[tgt_idx+i]);
            if(!isnan(temp))
            {
                dist += mask[i][j]*temp;
                den += mask[i][j];
            }
        }
        src_idx += imgM;
        tgt_idx += imgM;
    }

    /* channel B */
    src_idx = sj*imgM+si+2*numels;
    tgt_idx = tj*imgM+ti+2*numels; /* shift to the next layer */
    for(j=0; j<patchsize; ++j)
    {
        for(i=0; i<patchsize; ++i)
        {
            temp = mask[i][j]*fabs(img[src_idx+i]-img[tgt_idx+i]);
            if(!isnan(temp))
            {
                dist += mask[i][j]*temp;
                den += mask[i][j];
            }
        }
        src_idx += imgM;
        tgt_idx += imgM;
    }
    dist /= den;
    if(den<=threshold) dist = static_cast<double>(LPOS_VAL);
    return dist;
}


__inline double blockmatchgrayi_poisson(DATA_TYPE si, DATA_TYPE sj, DATA_TYPE ti, DATA_TYPE tj, double **mask, double *img, DATA_TYPE patchsize, DATA_TYPE imgM, DATA_TYPE imgN, double threshold)
{
    DATA_TYPE i, j, src_idx = sj*imgM+si, tgt_idx = tj*imgM+ti;
    double dist = 0.0, den = 0, temp;
    for(j=0; j<patchsize; ++j)
    {
        for(i=0; i<patchsize; ++i)
        {
            temp = (lgamma(img[src_idx+i]+img[tgt_idx+i]+0.5)/sqrt(lgamma(2*img[src_idx+i]+0.5)*lgamma(2*img[tgt_idx+i]+0.5)));
            if(!isnan(temp))
            {
                dist += mask[i][j]*temp;
                den += mask[i][j];
            }
            //mexPrintf("(si, sj) = (%d, %d), (ti, tj) = (%d, %d),|| i = %d, j= %d, || sval = %f, tval = %f, dist = %f\n", si,sj,ti,tj, i,j,img[src_idx+i],img[tgt_idx+i], (img[src_idx+i]-img[tgt_idx+i])*(img[src_idx+i]-img[tgt_idx+i]));
        }
        src_idx += imgM;
        tgt_idx += imgM;
    }
    dist /= den;
    if(den<=threshold) dist = static_cast<double>(LPOS_VAL);
    return dist;
}

__inline double blockmatchcolori_poisson(DATA_TYPE si, DATA_TYPE sj, DATA_TYPE ti, DATA_TYPE tj, double **mask, double *img, DATA_TYPE patchsize, DATA_TYPE imgM, DATA_TYPE imgN, double threshold)
{
    DATA_TYPE i, j, src_idx = sj*imgM+si, tgt_idx = tj*imgM+ti, numels = imgM*imgN;
    double dist = 0.0, den = 0, temp;
    /* channel R */
    for(j=0; j<patchsize; ++j)
    {
        for(i=0; i<patchsize; ++i)
        {
            temp = mask[i][j]*(lgamma(img[src_idx+i]+img[tgt_idx+i]+0.5)/sqrt(lgamma(2*img[src_idx+i]+0.5)*lgamma(2*img[tgt_idx+i]+0.5)));
            if(!isnan(temp))
            {
                dist += mask[i][j]*temp;
                den += mask[i][j];
            }
        }
        src_idx += imgM;
        tgt_idx += imgM;
    }

    /* channel G */
    src_idx = sj*imgM+si+numels;
    tgt_idx = tj*imgM+ti+numels; /* shift to the next layer */
    for(j=0; j<patchsize; ++j)
    {
        for(i=0; i<patchsize; ++i)
        {
            temp = mask[i][j]*(lgamma(img[src_idx+i]+img[tgt_idx+i]+0.5)/sqrt(lgamma(2*img[src_idx+i]+0.5)*lgamma(2*img[tgt_idx+i]+0.5)));
            if(!isnan(temp))
            {
                dist += mask[i][j]*temp;
                den += mask[i][j];
            }
        }
        src_idx += imgM;
        tgt_idx += imgM;
    }

    /* channel B */
    src_idx = sj*imgM+si+2*numels;
    tgt_idx = tj*imgM+ti+2*numels; /* shift to the next layer */
    for(j=0; j<patchsize; ++j)
    {
        for(i=0; i<patchsize; ++i)
        {
            temp = mask[i][j]*(lgamma(img[src_idx+i]+img[tgt_idx+i]+0.5)/sqrt(lgamma(2*img[src_idx+i]+0.5)*lgamma(2*img[tgt_idx+i]+0.5)));
            if(!isnan(temp))
            {
                dist += mask[i][j]*temp;
                den += mask[i][j];
            }
        }
        src_idx += imgM;
        tgt_idx += imgM;
    }
    dist /= den;
    if(den<=threshold) dist = static_cast<double>(LPOS_VAL);
    return dist;
}

