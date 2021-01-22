#ifndef __CONSTANTS_UTIL_H
#define __CONSTANTS_UTIL_H

// Physics constants
double gam, bc_mach, bc_alpha, bc_p, bc_r, bc_u, bc_e, dt;
double vortex_x0 = 5.0;
double vortex_y0 = 0.0;
double vortex_beta = 5.0;
double velX = 1.0;
double velY = 1.0;

// Utils
double ones[15] = {
  1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
};

int FMASK[15] = {0, 1, 2, 3, 4, 4, 8, 11, 13, 14, 14, 12, 9, 5, 0};

#endif
