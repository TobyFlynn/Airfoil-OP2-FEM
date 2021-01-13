inline void euler_fluxes(const double *q, double *F, double *G) {
  for(int i = 0; i < NUM_SOLUTION_PTS; i++) {
    double rho = q[i * 4 + 0];
    double rhou = q[i * 4 + 1];
    double rhov = q[i * 4 + 2];
    double ener = q[i * 4 + 3];
    // Calculate primative variables
    double u = rhou / rho;
    double v = rhov / rho;
    double p = (gam - 1) * (ener - 0.5 * (rhou * u + rhov * v));
    // Calculate flux
    F[i * 4 + 0] = rhou;
    F[i * 4 + 1] = rhou * u + p;
    F[i * 4 + 2] = rhov * u;
    F[i * 4 + 3] = u * (ener + p);
    G[i * 4 + 0] = rhov;
    G[i * 4 + 1] = rhou * v;
    G[i * 4 + 2] = rhov *v + p;
    G[i * 4 + 3] = v * (ener + p);
  }
}
