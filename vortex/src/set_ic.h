#include <cmath>

inline void set_ic(const double *t, const double *x, const double *y, double *q,
                   double *workingQ) {
  for(int i = 0; i < NUM_SOLUTION_PTS; i++) {
    double r = (x[i] - *t - vortex_x0) * (x[i] - *t - vortex_x0) + (y[i] - vortex_y0) * (y[i] - vortex_y0);
    r = sqrt(r);
    double u = 1.0 - vortex_beta * exp(1.0 - r * r) * ((y[i] - vortex_y0) / (2 * M_PI));
    double v = vortex_beta * exp(1.0 - r * r) * ((x[i] - vortex_x0) / (2 * M_PI));
    double rho = 1.0 - ((gam - 1) / (16 * gam * M_PI * M_PI)) * vortex_beta * vortex_beta * exp(2 * (1.0 - r * r));
    rho = pow(rho, 1.0 / (gam - 1));
    double p = pow(rho, gam);
    // u = 0.0;
    // v = 0.0;
    // rho = 1.0;
    // p = pow(rho, gam);
    q[i * 4]     = rho;
    q[i * 4 + 1] = rho * u;
    q[i * 4 + 2] = rho * v;
    // q[i * 4 + 3] = (p / (1.0 - gam)) + 0.5 * (rho * u + rho * v);
    q[i * 4 + 3] = (p / (gam - 1)) + 0.5 * rho * (u * u + v * v);
    workingQ[i * 4]     = q[i * 4];
    workingQ[i * 4 + 1] = q[i * 4 + 1];
    workingQ[i * 4 + 2] = q[i * 4 + 2];
    workingQ[i * 4 + 3] = q[i * 4 + 3];
    // cout << q[i * 4] << " " << q[i * 4 + 1] << " " << q[i * 4 + 2] << " " << q[i * 4 + 3] << endl;
  }
}
