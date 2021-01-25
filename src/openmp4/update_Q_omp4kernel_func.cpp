//
// auto-generated by op2.py
//

void update_Q_omp4_kernel(
  double *arg0,
  double *data1,
  int dat1size,
  double *data2,
  int dat2size,
  double *data3,
  int dat3size,
  double *data4,
  int dat4size,
  double *data5,
  int dat5size,
  int count,
  int num_teams,
  int nthread){

  double arg0_l = *arg0;
  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data1[0:dat1size],data2[0:dat2size],data3[0:dat3size],data4[0:dat4size],data5[0:dat5size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int n_op=0; n_op<count; n_op++ ){
    //variable mapping
    const double *dt = &arg0_l;
    double *q = &data1[60*n_op];
    const double *rk1 = &data2[60*n_op];
    const double *rk2 = &data3[60*n_op];
    const double *rk3 = &data4[60*n_op];
    double *workingQ = &data5[60*n_op];

    //inline function
    




    for(int i = 0; i < 4 * 15; i++) {
      q[i] = q[i] + (*dt) * (rk1[i]/ 6.0 + rk2[i] / 6.0 + 2.0 * rk3[i] / 3.0);
      workingQ[i] = q[i];
    }
    //end inline func
  }

  *arg0 = arg0_l;
}
