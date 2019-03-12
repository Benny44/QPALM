#include "mex.h"
#include <string.h>
#include "qpalm.h"
#include "cs.h"
#include "constants.h"

//Modes of operation
#define MODE_INIT "new"
#define MODE_DEFAULT_SETTINGS "default_settings"
#define MODE_RETRIEVE_CONSTANT "constant"
#define MODE_SETUP "setup"
#define MODE_SOLVE "solve"
#define MODE_DELETE "delete"

// QPALM work identifier
static QPALMWorkspace* qpalm_work = NULL;  

// all of the QPALM_INFO fieldnames as strings
const char* QPALM_INFO_FIELDS[] = {"iter",          //c_int
                                  "iter_out",       //c_int
                                  "status" ,        //char*
                                  "status_val" ,    //c_int
                                  "obj_val",        //c_float
                                  "pri_res_norm",   //c_float
                                  "dua_res_norm",   //c_float
                                  "dua2_res_norm",  //c_float
                                  "setup_time",     //c_float, only used if PROFILING
                                  "solve_time",     //c_float, only used if PROFILING
                                  "run_time"};      //c_float, only used if PROFILING


const char* QPALM_SETTINGS_FIELDS[] = {"max_iter",      //c_int
                                      "eps_abs",        //c_float
                                      "eps_rel",        //c_float
                                      "eps_abs_in",     //c_float
                                      "eps_rel_in",     //c_float
                                      "rho",            //c_float
                                      "eps_prim_inf",   //c_float
                                      "eps_dual_inf",   //c_float
                                      "theta",          //c_float
                                      "delta",          //c_float
                                      "tau_init",       //c_float
                                      "memory",         //c_int
                                      "proximal",       //c_int
                                      "gamma",          //c_int
                                      "gamma_upd",      //c_int
                                      "gamma_max",      //c_int
                                      "scaling",        //c_int
                                      "warm_start",     //c_int
                                      "verbose"};       //c_int

const char* CSC_FIELDS[] = {"nzmax",    //c_int
                            "m",        //c_int
                            "n",        //c_int
                            "p",        //c_int*
                            "i",        //c_int*
                            "x",        //c_float*
                            "nz"};      //c_int

const char* QPALM_DATA_FIELDS[] = {"n",     //c_int
                                   "m",     //c_int
                                   "Q",     //csc
                                   "A",     //csc
                                   "q",     //c_float*
                                   "bmin",  //c_float*
                                   "bmax"}; //c_float*

const char* QPALM_SCALING_FIELDS[] = {"D",       //c_float*
                                      "E",       //c_float*
                                      "Dinv",    //c_float*
                                      "Einv"};   //c_float*


const char* QPALM_WORKSPACE_FIELDS[] = {"data",
                                       "lbfgs",
                                       "scaling",
                                       "settings"};


// internal utility functions
void      castToDoubleArr(c_float *arr, double* arr_out, c_int len);
void      setToNaN(double* arr_out, c_int len);
void      copyMxStructToSettings(const mxArray*, QPALMSettings*);
mxArray*  copySettingsToMxStruct(QPALMSettings* settings);
mxArray*  copyInfoToMxStruct(QPALMInfo* info);
c_int*    copyToCintVector(mwIndex * vecData, c_int numel);
c_float*  copyToCfloatVector(double * vecData, c_int numel);



/**
 * The gateway function to QPALM
 *
 * qpalm_mex('default_settings');
 * qpalm_mex('constant', constant_name)
 * qpalm_mex('setup',n,m,Q,q,A,bmin,bmax,theSettings);
 * [out.x, out.y, out.prim_inf_cert, out.dual_inf_cert, out.info] = qpalm_mex('solve');
 * qpalm_mex('delete');
 *
 * @param nlhs Number of output arguments
 * @param plhs Array of output argument pointers
 * @param nrhs Number of input arguments
 * @param prhs Array of input argument pointers
 */
void mexFunction(int nlhs, mxArray * plhs [], int nrhs, const mxArray * prhs []) {

    // Get the command string
    char cmd[64];

    if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
		mexErrMsgTxt("First input should be a command string less than 64 characters long.");
    
    
    // report the default settings
    if (strcmp(cmd, MODE_DEFAULT_SETTINGS) == 0) {
        // Warn if other commands were ignored
        if (nrhs > 2)
            mexWarnMsgTxt("Default settings: unexpected number of arguments.");

        //Create a Settings structure in default form and report the results
        //Useful for external solver packages (e.g. Yalmip) that want to
        //know which solver settings are supported
        QPALMSettings* defaults = (QPALMSettings *)mxCalloc(1,sizeof(QPALMSettings));
        qpalm_set_default_settings(defaults);
        plhs[0] = copySettingsToMxStruct(defaults);
        mxFree(defaults);
        return;

    } else if (strcmp(cmd, MODE_DELETE) == 0) {
        
        //clean up the problem workspace
        if(qpalm_work != NULL){
            qpalm_cleanup(qpalm_work);
            qpalm_work = NULL;
        }
        // Warn if other commands were ignored
        if (nlhs != 0 || nrhs != 1)
            mexWarnMsgTxt("Delete: Unexpected arguments ignored.");
        return;


    } else if (strcmp(cmd, MODE_SETUP) == 0) {

        //throw an error if this is called more than once
        if(qpalm_work != NULL){
          mexErrMsgTxt("Solver is already initialized with problem data.");
        }

        //Create data and settings containers
        QPALMSettings* settings = (QPALMSettings *)mxCalloc(1,sizeof(QPALMSettings));
        QPALMData*     data     = (QPALMData *)mxCalloc(1,sizeof(QPALMData));

        // handle the problem data first.  Matlab-side
        // class wrapper is responsible for ensuring that
        // Q and A are sparse matrices,  everything
        // else is a dense vector and all inputs are
        // of compatible dimension

        // Get mxArrays
        const mxArray* Q  = prhs[3];
        const mxArray* q  = prhs[4];
        const mxArray* A  = prhs[5];
        const mxArray* bmin = prhs[6];
        const mxArray* bmax = prhs[7];

        // Create Data Structure
        data->n = (c_int) mxGetScalar(prhs[1]);
        data->m = (c_int) mxGetScalar(prhs[2]);
        data->q = copyToCfloatVector(mxGetPr(q), data->n);
        data->bmin = copyToCfloatVector(mxGetPr(bmin), data->m);
        data->bmax = copyToCfloatVector(mxGetPr(bmax), data->m);

        // Matrix Q:  nnz = Q->p[n]
        c_int * Qp = copyToCintVector(mxGetJc(Q), data->n + 1);
        c_int * Qi = copyToCintVector(mxGetIr(Q), Qp[data->n]);
        c_float * Qx = copyToCfloatVector(mxGetPr(Q), Qp[data->n]);
        data->Q  = csc_matrix(data->n, data->n, Qp[data->n], Qx, Qi, Qp);

        // Matrix A: nnz = A->p[n]
        c_int * Ap = copyToCintVector(mxGetJc(A), data->n + 1);
        c_int * Ai = copyToCintVector(mxGetIr(A), Ap[data->n]);
        c_float * Ax = copyToCfloatVector(mxGetPr(A), Ap[data->n]);
        data->A  = csc_matrix(data->m, data->n, Ap[data->n], Ax, Ai, Ap);

        // Create Settings
        const mxArray* mxSettings = prhs[8];
        if(mxIsEmpty(mxSettings)){
          // use defaults
          qpalm_set_default_settings(settings);
        } else {
          //populate settings structure from mxArray input
          copyMxStructToSettings(mxSettings, settings);
        }

        // Setup workspace

        qpalm_work = qpalm_setup(data, settings);
        if(qpalm_work == NULL){
           mexErrMsgTxt("Invalid problem setup");
         }

         //cleanup temporary structures
         // Data
         if (data->q) mxFree(data->q);
         if (data->bmin) mxFree(data->bmin);
         if (data->bmax) mxFree(data->bmax);
         if (Qx) mxFree(Qx);
         if (Qi) mxFree(Qi);
         if (Qp) mxFree(Qp);
         if (Ax) mxFree(Ax);
         if (Ai) mxFree(Ai);
         if (Ap) mxFree(Ap);
         mxFree(data);
         // Settings
         mxFree(settings);
        return;


    } else if (strcmp(cmd, MODE_SOLVE) == 0) { // SOLVE

        if (nlhs != 5 || nrhs != 1){
          mexErrMsgTxt("Solve : wrong number of inputs / outputs");
        }
        if(!qpalm_work){
            mexErrMsgTxt("No problem data has been given.");
        }
        // solve the problem
        qpalm_solve(qpalm_work);


        // Allocate space for solution
        // primal variables
        plhs[0] = mxCreateDoubleMatrix(qpalm_work->data->n,1,mxREAL);
        // dual variables
        plhs[1] = mxCreateDoubleMatrix(qpalm_work->data->m,1,mxREAL);
        // primal infeasibility certificate
        plhs[2] = mxCreateDoubleMatrix(qpalm_work->data->m,1,mxREAL);
        // dual infeasibility certificate
        plhs[3] = mxCreateDoubleMatrix(qpalm_work->data->n,1,mxREAL);

        //copy results to mxArray outputs
        //assume that five outputs will always
        //be returned to matlab-side class wrapper
        if ((qpalm_work->info->status_val != QPALM_PRIMAL_INFEASIBLE) &&
            (qpalm_work->info->status_val != QPALM_DUAL_INFEASIBLE)){

            //primal and dual solutions
            castToDoubleArr(qpalm_work->solution->x, mxGetPr(plhs[0]), qpalm_work->data->n);
            castToDoubleArr(qpalm_work->solution->y, mxGetPr(plhs[1]), qpalm_work->data->m);

            //infeasibility certificates -> NaN values
            setToNaN(mxGetPr(plhs[2]), qpalm_work->data->m);
            setToNaN(mxGetPr(plhs[3]), qpalm_work->data->n);

        } else if (qpalm_work->info->status_val == QPALM_PRIMAL_INFEASIBLE){ //primal infeasible

            //primal and dual solutions -> NaN values
            setToNaN(mxGetPr(plhs[0]), qpalm_work->data->n);
            setToNaN(mxGetPr(plhs[1]), qpalm_work->data->m);

            //primal infeasibility certificates
            castToDoubleArr(qpalm_work->delta_y, mxGetPr(plhs[2]), qpalm_work->data->m);

            //dual infeasibility certificates -> NaN values
            setToNaN(mxGetPr(plhs[3]), qpalm_work->data->n);

            // Set objective value to infinity
            qpalm_work->info->obj_val = mxGetInf();

        } else if (qpalm_work->info->status_val == QPALM_DUAL_INFEASIBLE) { //dual infeasible

            //primal and dual solutions -> NaN values
            setToNaN(mxGetPr(plhs[0]), qpalm_work->data->n);
            setToNaN(mxGetPr(plhs[1]), qpalm_work->data->m);

            //primal infeasibility certificates -> NaN values
            setToNaN(mxGetPr(plhs[2]), qpalm_work->data->m);

            //dual infeasibility certificates
            castToDoubleArr(qpalm_work->delta_x, mxGetPr(plhs[3]), qpalm_work->data->n);

            // Set objective value to -infinity
            qpalm_work->info->obj_val = -mxGetInf();
        }

        // if (qpalm_work->info->status_val == QPALM_NON_CVX) {
        //     qpalm_work->info->obj_val = mxGetNaN();
        // }

        plhs[4] = copyInfoToMxStruct(qpalm_work->info); // Info structure

        return;



    } else if (strcmp(cmd, MODE_RETRIEVE_CONSTANT) == 0) { // Return solver constants
         
        char constant[32];
        mxGetString(prhs[1], constant, sizeof(constant));

        if (!strcmp("QPALM_INFTY", constant)){
            plhs[0] = mxCreateDoubleScalar(QPALM_INFTY);
            return;
        }
        if (!strcmp("QPALM_NAN", constant)){
            plhs[0] = mxCreateDoubleScalar(mxGetNaN());
            return;
        }

        if (!strcmp("QPALM_SOLVED", constant)){
            plhs[0] = mxCreateDoubleScalar(QPALM_SOLVED);
            return;
        }

        if (!strcmp("QPALM_UNSOLVED", constant)){
            plhs[0] = mxCreateDoubleScalar(QPALM_UNSOLVED);
            return;
        }

        if (!strcmp("QPALM_PRIMAL_INFEASIBLE", constant)){
            plhs[0] = mxCreateDoubleScalar(QPALM_PRIMAL_INFEASIBLE);
            return;
        }

        if (!strcmp("QPALM_DUAL_INFEASIBLE", constant)){
            plhs[0] = mxCreateDoubleScalar(QPALM_DUAL_INFEASIBLE);
            return;
        }

        if (!strcmp("QPALM_MAX_ITER_REACHED", constant)){
            plhs[0] = mxCreateDoubleScalar(QPALM_MAX_ITER_REACHED);
            return;
        }

        if (!strcmp("QPALM_NON_CVX", constant)){
            plhs[0] = mxCreateDoubleScalar(QPALM_NON_CVX);
            return;
        }

        mexErrMsgTxt("Constant not recognized.");

        return;


    } else {
        //mexErrMsgIdAndTxt(PANOC_ERROR_RHS_OUT_OF_BOUNDS, "Invalid PANOC mode");
    }

}

void castToDoubleArr(c_float *arr, double* arr_out, c_int len) {
    for (c_int i = 0; i < len; i++) {
        arr_out[i] = (double)arr[i];
    }
}

void setToNaN(double* arr_out, c_int len){
    c_int i;
    for (i = 0; i < len; i++) {
        arr_out[i] = mxGetNaN();
    }
}

c_float*  copyToCfloatVector(double* vecData, c_int numel){
  // This memory needs to be freed!
  c_float* out = (c_float*)c_malloc(numel * sizeof(c_float));

  //copy data
  for(c_int i=0; i < numel; i++){
      out[i] = (c_float)vecData[i];
  }
  return out;

}

c_int* copyToCintVector(mwIndex* vecData, c_int numel){
  // This memory needs to be freed!
  c_int* out = (c_int*)c_malloc(numel * sizeof(c_int));

  //copy data
  for(c_int i=0; i < numel; i++){
      out[i] = (c_int)vecData[i];
  }
  return out;

}

mxArray* copyInfoToMxStruct(QPALMInfo* info){

  //create mxArray with the right number of fields
  int nfields  = sizeof(QPALM_INFO_FIELDS) / sizeof(QPALM_INFO_FIELDS[0]);
  mxArray* mxPtr = mxCreateStructMatrix(1,1,nfields,QPALM_INFO_FIELDS);

  //map the QPALM_INFO fields one at a time into mxArrays
  //matlab all numeric values as doubles
  mxSetField(mxPtr, 0, "iter",          mxCreateDoubleScalar(info->iter));
  mxSetField(mxPtr, 0, "iter_out",          mxCreateDoubleScalar(info->iter_out));
  mxSetField(mxPtr, 0, "status",        mxCreateString(info->status));
  mxSetField(mxPtr, 0, "status_val",    mxCreateDoubleScalar(info->status_val));
  mxSetField(mxPtr, 0, "obj_val",       mxCreateDoubleScalar(info->obj_val));
  mxSetField(mxPtr, 0, "pri_res_norm",       mxCreateDoubleScalar(info->pri_res_norm));
  mxSetField(mxPtr, 0, "dua_res_norm",       mxCreateDoubleScalar(info->dua_res_norm));
  mxSetField(mxPtr, 0, "dua2_res_norm",       mxCreateDoubleScalar(info->dua2_res_norm));

  #ifdef PROFILING
  //if not profiling, these fields will be empty
  mxSetField(mxPtr, 0, "setup_time",  mxCreateDoubleScalar(info->setup_time));
  mxSetField(mxPtr, 0, "solve_time",  mxCreateDoubleScalar(info->solve_time));
  mxSetField(mxPtr, 0, "run_time",    mxCreateDoubleScalar(info->run_time));
  #endif

  return mxPtr;

}

mxArray* copySettingsToMxStruct(QPALMSettings* settings){

  int nfields  = sizeof(QPALM_SETTINGS_FIELDS) / sizeof(QPALM_SETTINGS_FIELDS[0]);
  mxArray* mxPtr = mxCreateStructMatrix(1,1,nfields,QPALM_SETTINGS_FIELDS);

  //map the QPALM_SETTINGS fields one at a time into mxArrays
  //matlab handles everything as a double
  mxSetField(mxPtr, 0, "max_iter",        mxCreateDoubleScalar(settings->max_iter));
  mxSetField(mxPtr, 0, "eps_abs",         mxCreateDoubleScalar(settings->eps_abs));
  mxSetField(mxPtr, 0, "eps_rel",         mxCreateDoubleScalar(settings->eps_rel));
  mxSetField(mxPtr, 0, "eps_abs_in",      mxCreateDoubleScalar(settings->eps_abs_in));
  mxSetField(mxPtr, 0, "eps_rel_in",      mxCreateDoubleScalar(settings->eps_rel_in));
  mxSetField(mxPtr, 0, "rho",             mxCreateDoubleScalar(settings->rho));
  mxSetField(mxPtr, 0, "eps_prim_inf",    mxCreateDoubleScalar(settings->eps_prim_inf));
  mxSetField(mxPtr, 0, "eps_dual_inf",    mxCreateDoubleScalar(settings->eps_dual_inf));
  mxSetField(mxPtr, 0, "theta",           mxCreateDoubleScalar(settings->theta));
  mxSetField(mxPtr, 0, "delta",           mxCreateDoubleScalar(settings->delta));
  mxSetField(mxPtr, 0, "tau_init",        mxCreateDoubleScalar(settings->tau_init));
  mxSetField(mxPtr, 0, "memory",          mxCreateDoubleScalar(settings->memory));
  mxSetField(mxPtr, 0, "proximal",        mxCreateDoubleScalar(settings->proximal));
  mxSetField(mxPtr, 0, "gamma",           mxCreateDoubleScalar(settings->gamma));
  mxSetField(mxPtr, 0, "gamma_upd",       mxCreateDoubleScalar(settings->gamma_upd));
  mxSetField(mxPtr, 0, "gamma_max",       mxCreateDoubleScalar(settings->gamma_max));
  mxSetField(mxPtr, 0, "scaling",         mxCreateDoubleScalar(settings->scaling));
  mxSetField(mxPtr, 0, "warm_start",      mxCreateDoubleScalar(settings->warm_start));
  mxSetField(mxPtr, 0, "verbose",         mxCreateDoubleScalar(settings->verbose));

  return mxPtr;
}


void copyMxStructToSettings(const mxArray* mxPtr, QPALMSettings* settings){

  //this function assumes that only a complete and validated structure
  //will be passed.  matlab mex-side function is responsible for checking
  //structure validity

  //map the QPALM_SETTINGS fields one at a time into mxArrays
  //matlab handles everything as a double
  settings->max_iter                  = (c_int)mxGetScalar(mxGetField(mxPtr, 0, "max_iter"));
  settings->eps_abs                   = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "eps_abs"));
  settings->eps_rel                   = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "eps_rel"));
  settings->eps_abs_in                = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "eps_abs_in"));
  settings->eps_rel_in                = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "eps_rel_in"));
  settings->rho                       = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "rho"));
  settings->eps_prim_inf              = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "eps_dual_inf"));
  settings->eps_dual_inf              = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "eps_dual_inf"));
  settings->theta                     = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "theta"));
  settings->delta                     = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "delta"));
  settings->tau_init                  = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "tau_init"));
  settings->memory                    = (c_int)mxGetScalar(mxGetField(mxPtr, 0, "memory"));
  settings->proximal                  = (c_int)mxGetScalar(mxGetField(mxPtr, 0, "proximal"));
  settings->gamma                     = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "gamma"));  
  settings->gamma_upd                 = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "gamma_upd"));
  settings->gamma_max                 = (c_float)mxGetScalar(mxGetField(mxPtr, 0, "gamma_max"));
  settings->scaling                   = (c_int)mxGetScalar(mxGetField(mxPtr, 0, "scaling"));
  settings->warm_start                = (c_int)mxGetScalar(mxGetField(mxPtr, 0, "warm_start"));
  settings->verbose                   = (c_int)mxGetScalar(mxGetField(mxPtr, 0, "verbose"));

}