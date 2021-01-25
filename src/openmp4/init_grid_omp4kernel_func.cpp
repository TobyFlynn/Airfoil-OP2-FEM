//
// auto-generated by op2.py
//

#include <cblas.h>

void init_grid_omp4_kernel(
  int *map0,
  int map0size,
  double *data3,
  int dat3size,
  double *data4,
  int dat4size,
  double *data5,
  int dat5size,
  double *data6,
  int dat6size,
  double *data7,
  int dat7size,
  double *data8,
  int dat8size,
  double *data9,
  int dat9size,
  double *data10,
  int dat10size,
  double *data11,
  int dat11size,
  double *data12,
  int dat12size,
  double *data13,
  int dat13size,
  double *data14,
  int dat14size,
  double *data15,
  int dat15size,
  double *data16,
  int dat16size,
  double *data17,
  int dat17size,
  double *data0,
  int dat0size,
  int *col_reord,
  int set_size1,
  int start,
  int end,
  int num_teams,
  int nthread){

  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data3[0:dat3size],data4[0:dat4size],data5[0:dat5size],data6[0:dat6size],data7[0:dat7size],data8[0:dat8size],data9[0:dat9size],data10[0:dat10size],data11[0:dat11size],data12[0:dat12size],data13[0:dat13size],data14[0:dat14size],data15[0:dat15size],data16[0:dat16size],data17[0:dat17size]) \
    map(to: ones_ompkernel[:15], r_ompkernel[:15], s_ompkernel[:15], Dr_ompkernel[:225], Ds_ompkernel[:225], FMASK_ompkernel[:15])\
    map(to:col_reord[0:set_size1],map0[0:map0size],data0[0:dat0size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int e=start; e<end; e++ ){
    int n_op = col_reord[e];
    int map0idx;
    int map1idx;
    int map2idx;
    map0idx = map0[n_op + set_size1 * 0];
    map1idx = map0[n_op + set_size1 * 1];
    map2idx = map0[n_op + set_size1 * 2];

    //variable mapping
    const double *n0 = &data0[2 * map0idx];
    const double *n1 = &data0[2 * map1idx];
    const double *n2 = &data0[2 * map2idx];
    double *nodeX = &data3[3*n_op];
    double *nodeY = &data4[3*n_op];
    double *x = &data5[15*n_op];
    double *y = &data6[15*n_op];
    double *xr = &data7[15*n_op];
    double *yr = &data8[15*n_op];
    double *xs = &data9[15*n_op];
    double *ys = &data10[15*n_op];
    double *rx = &data11[15*n_op];
    double *ry = &data12[15*n_op];
    double *sx = &data13[15*n_op];
    double *sy = &data14[15*n_op];
    double *nx = &data15[15*n_op];
    double *ny = &data16[15*n_op];
    double *fscale = &data17[15*n_op];

    //inline function
    


    nodeX[0] = n0[0];
    nodeX[1] = n1[0];
    nodeX[2] = n2[0];
    nodeY[0] = n0[1];
    nodeY[1] = n1[1];
    nodeY[2] = n2[1];



    cblas_dcopy(15, ones_ompkernel, 1, x, 1);
    cblas_daxpy(15, 1.0, r_ompkernel, 1, x, 1);
    cblas_dscal(15, 0.5 * n1[0], x, 1);

    double temp[15];
    cblas_dcopy(15, ones_ompkernel, 1, temp, 1);
    cblas_daxpy(15, 1.0, s_ompkernel, 1, temp, 1);

    cblas_daxpy(15, 0.5 * n2[0], temp, 1, x, 1);

    cblas_dcopy(15, s_ompkernel, 1, temp, 1);
    cblas_daxpy(15, 1.0, r_ompkernel, 1, temp, 1);

    cblas_daxpy(15, -0.5 * n0[0], temp, 1, x, 1);



    cblas_dcopy(15, ones_ompkernel, 1, y, 1);
    cblas_daxpy(15, 1.0, r_ompkernel, 1, y, 1);
    cblas_dscal(15, 0.5 * n1[1], y, 1);

    cblas_dcopy(15, ones_ompkernel, 1, temp, 1);
    cblas_daxpy(15, 1.0, s_ompkernel, 1, temp, 1);

    cblas_daxpy(15, 0.5 * n2[1], temp, 1, y, 1);

    cblas_dcopy(15, s_ompkernel, 1, temp, 1);
    cblas_daxpy(15, 1.0, r_ompkernel, 1, temp, 1);

    cblas_daxpy(15, -0.5 * n0[1], temp, 1, y, 1);


    cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, 1.0, Dr_ompkernel, 15, x, 1, 0.0, xr, 1);

    cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, 1.0, Ds_ompkernel, 15, x, 1, 0.0, xs, 1);

    cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, 1.0, Dr_ompkernel, 15, y, 1, 0.0, yr, 1);

    cblas_dgemv(CblasRowMajor, CblasNoTrans, 15, 15, 1.0, Ds_ompkernel, 15, y, 1, 0.0, ys, 1);

    double J[15];
    for(int i = 0; i < 15; i++) {
      J[i] = -xs[i] * yr[i] + xr[i] * ys[i];
    }

    for(int i = 0; i < 15; i++) {
      rx[i] = ys[i] / J[i];
      sx[i] = -yr[i] / J[i];
      ry[i] = -xs[i] / J[i];
      sy[i] = xr[i] / J[i];
    }


    for(int i = 0; i < 5; i++) {
      nx[i] = yr[FMASK_ompkernel[i]];
      ny[i] = -xr[FMASK_ompkernel[i]];
    }

    for(int i = 0; i < 5; i++) {
      nx[5 + i] = ys[FMASK_ompkernel[5 + i]] - yr[FMASK_ompkernel[5 + i]];
      ny[5 + i] = xr[FMASK_ompkernel[5 + i]] - xs[FMASK_ompkernel[5 + i]];
    }

    for(int i = 0; i < 5; i++) {
      nx[2 * 5 + i] = -ys[FMASK_ompkernel[2 * 5 + i]];
      ny[2 * 5 + i] = xs[FMASK_ompkernel[2 * 5 + i]];
    }

    for(int i = 0; i < 3 * 5; i++) {
      double sJ = sqrt(nx[i] * nx[i] + ny[i] * ny[i]);
      nx[i] = nx[i] / sJ;
      ny[i] = ny[i] / sJ;
      fscale[i] = sJ / J[FMASK_ompkernel[i]];
    }
    //end inline func
  }

}
