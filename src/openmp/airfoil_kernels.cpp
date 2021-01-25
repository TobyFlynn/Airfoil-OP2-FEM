//
// auto-generated by op2.py
//

#ifdef _OPENMP
  #include <omp.h>
#endif

// global constants
extern double gam;
extern double bc_mach;
extern double bc_alpha;
extern double bc_p;
extern double bc_r;
extern double bc_u;
extern double bc_v;
extern double bc_e;
extern double ones[15];
extern double r[15];
extern double s[15];
extern double Dr[225];
extern double Ds[225];
extern double Drw[225];
extern double Dsw[225];
extern int FMASK[15];
extern double LIFT[225];

// header
#include "op_lib_cpp.h"

// user kernel files
#include "init_grid_kernel.cpp"
#include "set_ic_kernel.cpp"
#include "flux_zero_kernel.cpp"
#include "calc_dt_kernel.cpp"
#include "get_neighbour_q_kernel.cpp"
#include "get_bedge_q_kernel.cpp"
#include "euler_rhs_kernel.cpp"
#include "set_workingQ_kernel.cpp"
#include "update_Q_kernel.cpp"
