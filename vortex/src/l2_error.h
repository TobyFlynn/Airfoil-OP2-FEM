#include <cmath>

inline void l2_error(const double *x, const double *y, const double *q,
                     const double *t, double *error_r, double *error_ru,
                     double *error_rv, double *error_ener, double *error) {
  *error_r    = 0.0;
  *error_ru   = 0.0;
  *error_rv   = 0.0;
  *error_ener = 0.0;
  for(int i = 0; i < NUM_SOLUTION_PTS; i++) {
    double r = (x[i] - *t - vortex_x0) * (x[i] - *t - vortex_x0) + (y[i] - vortex_y0) * (y[i] - vortex_y0);
    r = sqrt(r);
    double u = 1.0 - vortex_beta * exp(1.0 - r * r) * ((y[i] - vortex_y0) / (2 * M_PI));
    double v = vortex_beta * exp(1.0 - r * r) * ((x[i] - vortex_x0) / (2 * M_PI));
    double rho = 1.0 - ((gam - 1) / (16 * gam * M_PI * M_PI)) * vortex_beta * vortex_beta * exp(2 * (1.0 - r * r));
    rho = pow(rho, 1.0 / (gam - 1));
    double p = pow(rho, gam);
    // double ener = (p / (1.0 - gam)) + 0.5 * (rho * u + rho * v);
    double ener = (p / (gam - 1)) + 0.5 * rho * (u * u + v * v);
    *error_r    += (q[i * 4] - rho) * (q[i * 4] - rho);
    *error_ru   += (q[i * 4 + 1] - rho * u) * (q[i * 4 + 1] - rho * u);
    *error_rv   += (q[i * 4 + 2] - rho * v) * (q[i * 4 + 2] - rho * v);
    *error_ener += (q[i * 4 + 3] - ener) * (q[i * 4 + 3] - ener);
    error[i * 4] = (q[i * 4] - rho) * (q[i * 4] - rho);
    error[i * 4 + 1] = (q[i * 4 + 1] - rho * u) * (q[i * 4 + 1] - rho * u);
    error[i * 4 + 2] = (q[i * 4 + 2] - rho * v) * (q[i * 4 + 2] - rho * v);
    error[i * 4 + 3] = (q[i * 4 + 3] - ener) * (q[i * 4 + 3] - ener);
  }
}
