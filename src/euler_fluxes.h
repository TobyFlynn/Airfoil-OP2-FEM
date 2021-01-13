inline void euler_fluxes(const double *q, double *F, double *G) {
  double rho = q[0];
  double rhou = q[1];
  double rhov = q[2];
  double ener = q[3];
  // Calculate primative variables
  double u = rhou / rho;
  double v = rhov / rho;
  double p = (gam - 1) * (ener - 0.5 * (rhou * u + rhov * v));
  // Calculate flux
  F[0] = rhou;
  F[1] = rhou * u + p;
  F[2] = rhov * u;
  F[3] = u * (ener + p);
  G[0] = rhov;
  G[1] = rhou * v;
  G[2] = rhov *v + p;
  G[3] = v * (ener + p);
}
