//
// auto-generated by op2.py
//

// global constants
double gam_ompkernel;
double bc_mach_ompkernel;
double bc_alpha_ompkernel;
double bc_p_ompkernel;
double bc_r_ompkernel;
double bc_u_ompkernel;
double bc_v_ompkernel;
double bc_e_ompkernel;
double ones_ompkernel[15];
double r_ompkernel[15];
double s_ompkernel[15];
double Dr_ompkernel[225];
double Ds_ompkernel[225];
double Drw_ompkernel[225];
double Dsw_ompkernel[225];
int FMASK_ompkernel[15];
double LIFT_ompkernel[225];

// header
#include "op_lib_cpp.h"

void op_decl_const_char(int dim, char const *type,
  int size, char *dat, char const *name){
  if(!strcmp(name, "gam")) {
    memcpy(&gam_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:gam_ompkernel)
  } else if(!strcmp(name, "bc_mach")) {
    memcpy(&bc_mach_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:bc_mach_ompkernel)
  } else if(!strcmp(name, "bc_alpha")) {
    memcpy(&bc_alpha_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:bc_alpha_ompkernel)
  } else if(!strcmp(name, "bc_p")) {
    memcpy(&bc_p_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:bc_p_ompkernel)
  } else if(!strcmp(name, "bc_r")) {
    memcpy(&bc_r_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:bc_r_ompkernel)
  } else if(!strcmp(name, "bc_u")) {
    memcpy(&bc_u_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:bc_u_ompkernel)
  } else if(!strcmp(name, "bc_v")) {
    memcpy(&bc_v_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:bc_v_ompkernel)
  } else if(!strcmp(name, "bc_e")) {
    memcpy(&bc_e_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:bc_e_ompkernel)
  } else if(!strcmp(name, "ones")) {
    memcpy(ones_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:ones_ompkernel[:15])
  } else if(!strcmp(name, "r")) {
    memcpy(r_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:r_ompkernel[:15])
  } else if(!strcmp(name, "s")) {
    memcpy(s_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:s_ompkernel[:15])
  } else if(!strcmp(name, "Dr")) {
    memcpy(Dr_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:Dr_ompkernel[:225])
  } else if(!strcmp(name, "Ds")) {
    memcpy(Ds_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:Ds_ompkernel[:225])
  } else if(!strcmp(name, "Drw")) {
    memcpy(Drw_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:Drw_ompkernel[:225])
  } else if(!strcmp(name, "Dsw")) {
    memcpy(Dsw_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:Dsw_ompkernel[:225])
  } else if(!strcmp(name, "FMASK")) {
    memcpy(FMASK_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:FMASK_ompkernel[:15])
  } else if(!strcmp(name, "LIFT")) {
    memcpy(LIFT_ompkernel, dat, dim*size);
  #pragma omp target enter data map(to:LIFT_ompkernel[:225])
  }
}
// user kernel files
#include "init_grid_omp4kernel_func.cpp"
#include "set_ic_omp4kernel_func.cpp"
#include "calc_dt_omp4kernel_func.cpp"
#include "get_neighbour_q_omp4kernel_func.cpp"
#include "get_bedge_q_omp4kernel_func.cpp"
#include "internal_fluxes_omp4kernel_func.cpp"
#include "euler_rhs_omp4kernel_func.cpp"
#include "set_workingQ_omp4kernel_func.cpp"
#include "update_Q_omp4kernel_func.cpp"
