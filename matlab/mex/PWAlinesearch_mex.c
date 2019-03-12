#include "mex.h"

#define IS_REAL_DENSE_VEC(P) ((mxGetNumberOfDimensions(P) == 1 || \
    (mxGetNumberOfDimensions(P) == 2 && (mxGetN(P) == 1 || mxGetM(P) == 1))) && \
    !mxIsSparse(P) && mxIsDouble(P))
#define IS_REAL_SCALAR(P) (IS_REAL_DENSE_VEC(P) && mxGetNumberOfElements(P) == 1)
#define IS_INT64_SCALAR(P) ((mxGetNumberOfDimensions(P) == 1 || \
    (mxGetNumberOfDimensions(P) == 2 && (mxGetN(P) == 1 && mxGetM(P) == 1))) && \
    mxGetNumberOfElements(P) == 1 && mxIsInt64(P))

/**
 * Array (with indices) to sort
 */
typedef struct array_element  {
  double x; ///< Array element
  int   i; ///< Index
} array_element;

int compare (const void * a, const void * b)
{
    double f = ((struct array_element*)a)->x;
    double s = ((struct array_element*)b)->x;
    if (f > s) return  1;
    if (f < s) return -1;
    return 0;
}

void mexFunction(int nlhs, mxArray * plhs [], int nrhs, const mxArray * prhs []) {

    if (nrhs != 5) {
        mexErrMsgTxt("PWAlinesearch: you should provide exactly 5 arguments.");
        return;
    }
    if (nlhs != 1) {
        mexErrMsgTxt("PWAlinesearch: you should expect only one output argument.");
        return;
    }
    if (!IS_REAL_SCALAR(prhs[0])) {
        mexErrMsgTxt("PWAlinesearch: 1st argument must be a double scalar.");
        return;
    }
    if (!IS_REAL_SCALAR(prhs[1])) {
        mexErrMsgTxt("PWAlinesearch: 2nd argument must be a double scalar.");
        return;
    }
    if (!IS_REAL_DENSE_VEC(prhs[2])) {
        mexErrMsgTxt("PWAlinesearch: 3rd argument must be a double, dense vector.");
        return;
    }
    if (!IS_REAL_DENSE_VEC(prhs[3])) {
        mexErrMsgTxt("PWAlinesearch: 4th argument must be a double, dense vector.");
        return;
    }
    if (!IS_INT64_SCALAR(prhs[4])) {
        mexErrMsgTxt("PWAlinesearch: 5th argument must be an int scalar.");
        return;
    }
    

    const double eta = mxGetScalar(prhs[0]);
    const double beta = mxGetScalar(prhs[1]);
    const double *delta = (double*)mxGetPr(prhs[2]);
    const double *alpha = (double*)mxGetPr(prhs[3]);
    int m = (int)mxGetScalar(prhs[4]);

    array_element *s =(array_element*) mxCalloc(m, sizeof(array_element));
    int *index_P = (int*) mxCalloc(m, sizeof(int));
    int *index_L = (int*) mxCalloc(m, sizeof(int));
    int *index_J = (int*) mxCalloc(m, sizeof(int));

    double ad; // alpha/delta
    int i = 0;
    int ns = 0;
    for (; i < m; i++) {
        ad = alpha[i]/delta[i];
        if (ad > 0) {
            s[ns].x = alpha[i]/delta[i];
            s[ns].i = i;
            ns++;
            index_L[i] = 1;
        } else {
            index_L[i] = 0;
        }
        index_P[i] = (delta[i] > 0);
        index_J[i] = ((index_P[i] + index_L[i]) == 1);
    }
        
    double a, b;
    a = eta;
    b = beta;
    for (i=0; i<m; i++) {
        a += index_J[i]*delta[i]*delta[i];
        b -= index_J[i]*delta[i]*alpha[i];
    }

    qsort(s, ns, sizeof(array_element), compare);

    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double *tau; 
    tau = mxGetPr(plhs[0]);

    if (ns == 0 || a*s[0].x+b > 0) {
        mxFree(s);
        mxFree(index_P);
        mxFree(index_L);
        mxFree(index_J);
        *tau = -b/a;
        return;
    }
    
    i = 0;
    int iz;

    while (i < ns-1) {
        iz = s[i].i;
        if (index_P[iz]) {
            a += delta[iz]*delta[iz];
            b -= delta[iz]*alpha[iz];
        } else {
            a -= delta[iz]*delta[iz];
            b += delta[iz]*alpha[iz];
        }
        i++;
        if (a*s[i].x+b > 0) {
            mxFree(s);
            mxFree(index_P);
            mxFree(index_L);
            mxFree(index_J);
            *tau = -b/a;
            return;
        }
    }
    iz = s[i].i;
    if (index_P[iz]) {
        a += delta[iz]*delta[iz];
        b -= delta[iz]*alpha[iz];
    }
    mxFree(s);
    mxFree(index_P);
    mxFree(index_L);
    mxFree(index_J);
    *tau = -b/a;
    return;

}