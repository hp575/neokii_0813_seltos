#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 5.991464547107981;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with sympy 1.8                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_5096038058135125524) {
   out_5096038058135125524[0] = delta_x[0] + nom_x[0];
   out_5096038058135125524[1] = delta_x[1] + nom_x[1];
   out_5096038058135125524[2] = delta_x[2] + nom_x[2];
   out_5096038058135125524[3] = delta_x[3] + nom_x[3];
   out_5096038058135125524[4] = delta_x[4] + nom_x[4];
   out_5096038058135125524[5] = delta_x[5] + nom_x[5];
   out_5096038058135125524[6] = delta_x[6] + nom_x[6];
   out_5096038058135125524[7] = delta_x[7] + nom_x[7];
   out_5096038058135125524[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_8373756370158955549) {
   out_8373756370158955549[0] = -nom_x[0] + true_x[0];
   out_8373756370158955549[1] = -nom_x[1] + true_x[1];
   out_8373756370158955549[2] = -nom_x[2] + true_x[2];
   out_8373756370158955549[3] = -nom_x[3] + true_x[3];
   out_8373756370158955549[4] = -nom_x[4] + true_x[4];
   out_8373756370158955549[5] = -nom_x[5] + true_x[5];
   out_8373756370158955549[6] = -nom_x[6] + true_x[6];
   out_8373756370158955549[7] = -nom_x[7] + true_x[7];
   out_8373756370158955549[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_7930719249271816004) {
   out_7930719249271816004[0] = 1.0;
   out_7930719249271816004[1] = 0.0;
   out_7930719249271816004[2] = 0.0;
   out_7930719249271816004[3] = 0.0;
   out_7930719249271816004[4] = 0.0;
   out_7930719249271816004[5] = 0.0;
   out_7930719249271816004[6] = 0.0;
   out_7930719249271816004[7] = 0.0;
   out_7930719249271816004[8] = 0.0;
   out_7930719249271816004[9] = 0.0;
   out_7930719249271816004[10] = 1.0;
   out_7930719249271816004[11] = 0.0;
   out_7930719249271816004[12] = 0.0;
   out_7930719249271816004[13] = 0.0;
   out_7930719249271816004[14] = 0.0;
   out_7930719249271816004[15] = 0.0;
   out_7930719249271816004[16] = 0.0;
   out_7930719249271816004[17] = 0.0;
   out_7930719249271816004[18] = 0.0;
   out_7930719249271816004[19] = 0.0;
   out_7930719249271816004[20] = 1.0;
   out_7930719249271816004[21] = 0.0;
   out_7930719249271816004[22] = 0.0;
   out_7930719249271816004[23] = 0.0;
   out_7930719249271816004[24] = 0.0;
   out_7930719249271816004[25] = 0.0;
   out_7930719249271816004[26] = 0.0;
   out_7930719249271816004[27] = 0.0;
   out_7930719249271816004[28] = 0.0;
   out_7930719249271816004[29] = 0.0;
   out_7930719249271816004[30] = 1.0;
   out_7930719249271816004[31] = 0.0;
   out_7930719249271816004[32] = 0.0;
   out_7930719249271816004[33] = 0.0;
   out_7930719249271816004[34] = 0.0;
   out_7930719249271816004[35] = 0.0;
   out_7930719249271816004[36] = 0.0;
   out_7930719249271816004[37] = 0.0;
   out_7930719249271816004[38] = 0.0;
   out_7930719249271816004[39] = 0.0;
   out_7930719249271816004[40] = 1.0;
   out_7930719249271816004[41] = 0.0;
   out_7930719249271816004[42] = 0.0;
   out_7930719249271816004[43] = 0.0;
   out_7930719249271816004[44] = 0.0;
   out_7930719249271816004[45] = 0.0;
   out_7930719249271816004[46] = 0.0;
   out_7930719249271816004[47] = 0.0;
   out_7930719249271816004[48] = 0.0;
   out_7930719249271816004[49] = 0.0;
   out_7930719249271816004[50] = 1.0;
   out_7930719249271816004[51] = 0.0;
   out_7930719249271816004[52] = 0.0;
   out_7930719249271816004[53] = 0.0;
   out_7930719249271816004[54] = 0.0;
   out_7930719249271816004[55] = 0.0;
   out_7930719249271816004[56] = 0.0;
   out_7930719249271816004[57] = 0.0;
   out_7930719249271816004[58] = 0.0;
   out_7930719249271816004[59] = 0.0;
   out_7930719249271816004[60] = 1.0;
   out_7930719249271816004[61] = 0.0;
   out_7930719249271816004[62] = 0.0;
   out_7930719249271816004[63] = 0.0;
   out_7930719249271816004[64] = 0.0;
   out_7930719249271816004[65] = 0.0;
   out_7930719249271816004[66] = 0.0;
   out_7930719249271816004[67] = 0.0;
   out_7930719249271816004[68] = 0.0;
   out_7930719249271816004[69] = 0.0;
   out_7930719249271816004[70] = 1.0;
   out_7930719249271816004[71] = 0.0;
   out_7930719249271816004[72] = 0.0;
   out_7930719249271816004[73] = 0.0;
   out_7930719249271816004[74] = 0.0;
   out_7930719249271816004[75] = 0.0;
   out_7930719249271816004[76] = 0.0;
   out_7930719249271816004[77] = 0.0;
   out_7930719249271816004[78] = 0.0;
   out_7930719249271816004[79] = 0.0;
   out_7930719249271816004[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_2019501809289876349) {
   out_2019501809289876349[0] = state[0];
   out_2019501809289876349[1] = state[1];
   out_2019501809289876349[2] = state[2];
   out_2019501809289876349[3] = state[3];
   out_2019501809289876349[4] = state[4];
   out_2019501809289876349[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_2019501809289876349[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_2019501809289876349[7] = state[7];
   out_2019501809289876349[8] = state[8];
}
void F_fun(double *state, double dt, double *out_3574153591255788280) {
   out_3574153591255788280[0] = 1;
   out_3574153591255788280[1] = 0;
   out_3574153591255788280[2] = 0;
   out_3574153591255788280[3] = 0;
   out_3574153591255788280[4] = 0;
   out_3574153591255788280[5] = 0;
   out_3574153591255788280[6] = 0;
   out_3574153591255788280[7] = 0;
   out_3574153591255788280[8] = 0;
   out_3574153591255788280[9] = 0;
   out_3574153591255788280[10] = 1;
   out_3574153591255788280[11] = 0;
   out_3574153591255788280[12] = 0;
   out_3574153591255788280[13] = 0;
   out_3574153591255788280[14] = 0;
   out_3574153591255788280[15] = 0;
   out_3574153591255788280[16] = 0;
   out_3574153591255788280[17] = 0;
   out_3574153591255788280[18] = 0;
   out_3574153591255788280[19] = 0;
   out_3574153591255788280[20] = 1;
   out_3574153591255788280[21] = 0;
   out_3574153591255788280[22] = 0;
   out_3574153591255788280[23] = 0;
   out_3574153591255788280[24] = 0;
   out_3574153591255788280[25] = 0;
   out_3574153591255788280[26] = 0;
   out_3574153591255788280[27] = 0;
   out_3574153591255788280[28] = 0;
   out_3574153591255788280[29] = 0;
   out_3574153591255788280[30] = 1;
   out_3574153591255788280[31] = 0;
   out_3574153591255788280[32] = 0;
   out_3574153591255788280[33] = 0;
   out_3574153591255788280[34] = 0;
   out_3574153591255788280[35] = 0;
   out_3574153591255788280[36] = 0;
   out_3574153591255788280[37] = 0;
   out_3574153591255788280[38] = 0;
   out_3574153591255788280[39] = 0;
   out_3574153591255788280[40] = 1;
   out_3574153591255788280[41] = 0;
   out_3574153591255788280[42] = 0;
   out_3574153591255788280[43] = 0;
   out_3574153591255788280[44] = 0;
   out_3574153591255788280[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_3574153591255788280[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_3574153591255788280[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_3574153591255788280[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_3574153591255788280[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_3574153591255788280[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_3574153591255788280[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_3574153591255788280[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_3574153591255788280[53] = -9.8000000000000007*dt;
   out_3574153591255788280[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_3574153591255788280[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_3574153591255788280[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3574153591255788280[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3574153591255788280[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_3574153591255788280[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_3574153591255788280[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_3574153591255788280[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3574153591255788280[62] = 0;
   out_3574153591255788280[63] = 0;
   out_3574153591255788280[64] = 0;
   out_3574153591255788280[65] = 0;
   out_3574153591255788280[66] = 0;
   out_3574153591255788280[67] = 0;
   out_3574153591255788280[68] = 0;
   out_3574153591255788280[69] = 0;
   out_3574153591255788280[70] = 1;
   out_3574153591255788280[71] = 0;
   out_3574153591255788280[72] = 0;
   out_3574153591255788280[73] = 0;
   out_3574153591255788280[74] = 0;
   out_3574153591255788280[75] = 0;
   out_3574153591255788280[76] = 0;
   out_3574153591255788280[77] = 0;
   out_3574153591255788280[78] = 0;
   out_3574153591255788280[79] = 0;
   out_3574153591255788280[80] = 1;
}
void h_25(double *state, double *unused, double *out_1634828461630621151) {
   out_1634828461630621151[0] = state[6];
}
void H_25(double *state, double *unused, double *out_2197806516232288857) {
   out_2197806516232288857[0] = 0;
   out_2197806516232288857[1] = 0;
   out_2197806516232288857[2] = 0;
   out_2197806516232288857[3] = 0;
   out_2197806516232288857[4] = 0;
   out_2197806516232288857[5] = 0;
   out_2197806516232288857[6] = 1;
   out_2197806516232288857[7] = 0;
   out_2197806516232288857[8] = 0;
}
void h_24(double *state, double *unused, double *out_8272163805170178838) {
   out_8272163805170178838[0] = state[4];
   out_8272163805170178838[1] = state[5];
}
void H_24(double *state, double *unused, double *out_25156917226789291) {
   out_25156917226789291[0] = 0;
   out_25156917226789291[1] = 0;
   out_25156917226789291[2] = 0;
   out_25156917226789291[3] = 0;
   out_25156917226789291[4] = 1;
   out_25156917226789291[5] = 0;
   out_25156917226789291[6] = 0;
   out_25156917226789291[7] = 0;
   out_25156917226789291[8] = 0;
   out_25156917226789291[9] = 0;
   out_25156917226789291[10] = 0;
   out_25156917226789291[11] = 0;
   out_25156917226789291[12] = 0;
   out_25156917226789291[13] = 0;
   out_25156917226789291[14] = 1;
   out_25156917226789291[15] = 0;
   out_25156917226789291[16] = 0;
   out_25156917226789291[17] = 0;
}
void h_30(double *state, double *unused, double *out_6049964250792746476) {
   out_6049964250792746476[0] = state[4];
}
void H_30(double *state, double *unused, double *out_4716139474739537484) {
   out_4716139474739537484[0] = 0;
   out_4716139474739537484[1] = 0;
   out_4716139474739537484[2] = 0;
   out_4716139474739537484[3] = 0;
   out_4716139474739537484[4] = 1;
   out_4716139474739537484[5] = 0;
   out_4716139474739537484[6] = 0;
   out_4716139474739537484[7] = 0;
   out_4716139474739537484[8] = 0;
}
void h_26(double *state, double *unused, double *out_4192496911818067350) {
   out_4192496911818067350[0] = state[7];
}
void H_26(double *state, double *unused, double *out_1543696802641767367) {
   out_1543696802641767367[0] = 0;
   out_1543696802641767367[1] = 0;
   out_1543696802641767367[2] = 0;
   out_1543696802641767367[3] = 0;
   out_1543696802641767367[4] = 0;
   out_1543696802641767367[5] = 0;
   out_1543696802641767367[6] = 0;
   out_1543696802641767367[7] = 1;
   out_1543696802641767367[8] = 0;
}
void h_27(double *state, double *unused, double *out_8015606147812933325) {
   out_8015606147812933325[0] = state[3];
}
void H_27(double *state, double *unused, double *out_4504653125695744252) {
   out_4504653125695744252[0] = 0;
   out_4504653125695744252[1] = 0;
   out_4504653125695744252[2] = 0;
   out_4504653125695744252[3] = 1;
   out_4504653125695744252[4] = 0;
   out_4504653125695744252[5] = 0;
   out_4504653125695744252[6] = 0;
   out_4504653125695744252[7] = 0;
   out_4504653125695744252[8] = 0;
}
void h_29(double *state, double *unused, double *out_8172792816164062470) {
   out_8172792816164062470[0] = state[1];
}
void H_29(double *state, double *unused, double *out_1819658469580927157) {
   out_1819658469580927157[0] = 0;
   out_1819658469580927157[1] = 1;
   out_1819658469580927157[2] = 0;
   out_1819658469580927157[3] = 0;
   out_1819658469580927157[4] = 0;
   out_1819658469580927157[5] = 0;
   out_1819658469580927157[6] = 0;
   out_1819658469580927157[7] = 0;
   out_1819658469580927157[8] = 0;
}
void h_28(double *state, double *unused, double *out_5938536161097710960) {
   out_5938536161097710960[0] = state[5];
   out_5938536161097710960[1] = state[6];
}
void H_28(double *state, double *unused, double *out_1495157348996847574) {
   out_1495157348996847574[0] = 0;
   out_1495157348996847574[1] = 0;
   out_1495157348996847574[2] = 0;
   out_1495157348996847574[3] = 0;
   out_1495157348996847574[4] = 0;
   out_1495157348996847574[5] = 1;
   out_1495157348996847574[6] = 0;
   out_1495157348996847574[7] = 0;
   out_1495157348996847574[8] = 0;
   out_1495157348996847574[9] = 0;
   out_1495157348996847574[10] = 0;
   out_1495157348996847574[11] = 0;
   out_1495157348996847574[12] = 0;
   out_1495157348996847574[13] = 0;
   out_1495157348996847574[14] = 0;
   out_1495157348996847574[15] = 1;
   out_1495157348996847574[16] = 0;
   out_1495157348996847574[17] = 0;
}
void h_31(double *state, double *unused, double *out_1792015129981750296) {
   out_1792015129981750296[0] = state[8];
}
void H_31(double *state, double *unused, double *out_2169904904875118843) {
   out_2169904904875118843[0] = 0;
   out_2169904904875118843[1] = 0;
   out_2169904904875118843[2] = 0;
   out_2169904904875118843[3] = 0;
   out_2169904904875118843[4] = 0;
   out_2169904904875118843[5] = 0;
   out_2169904904875118843[6] = 0;
   out_2169904904875118843[7] = 0;
   out_2169904904875118843[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_5096038058135125524) {
  err_fun(nom_x, delta_x, out_5096038058135125524);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_8373756370158955549) {
  inv_err_fun(nom_x, true_x, out_8373756370158955549);
}
void car_H_mod_fun(double *state, double *out_7930719249271816004) {
  H_mod_fun(state, out_7930719249271816004);
}
void car_f_fun(double *state, double dt, double *out_2019501809289876349) {
  f_fun(state,  dt, out_2019501809289876349);
}
void car_F_fun(double *state, double dt, double *out_3574153591255788280) {
  F_fun(state,  dt, out_3574153591255788280);
}
void car_h_25(double *state, double *unused, double *out_1634828461630621151) {
  h_25(state, unused, out_1634828461630621151);
}
void car_H_25(double *state, double *unused, double *out_2197806516232288857) {
  H_25(state, unused, out_2197806516232288857);
}
void car_h_24(double *state, double *unused, double *out_8272163805170178838) {
  h_24(state, unused, out_8272163805170178838);
}
void car_H_24(double *state, double *unused, double *out_25156917226789291) {
  H_24(state, unused, out_25156917226789291);
}
void car_h_30(double *state, double *unused, double *out_6049964250792746476) {
  h_30(state, unused, out_6049964250792746476);
}
void car_H_30(double *state, double *unused, double *out_4716139474739537484) {
  H_30(state, unused, out_4716139474739537484);
}
void car_h_26(double *state, double *unused, double *out_4192496911818067350) {
  h_26(state, unused, out_4192496911818067350);
}
void car_H_26(double *state, double *unused, double *out_1543696802641767367) {
  H_26(state, unused, out_1543696802641767367);
}
void car_h_27(double *state, double *unused, double *out_8015606147812933325) {
  h_27(state, unused, out_8015606147812933325);
}
void car_H_27(double *state, double *unused, double *out_4504653125695744252) {
  H_27(state, unused, out_4504653125695744252);
}
void car_h_29(double *state, double *unused, double *out_8172792816164062470) {
  h_29(state, unused, out_8172792816164062470);
}
void car_H_29(double *state, double *unused, double *out_1819658469580927157) {
  H_29(state, unused, out_1819658469580927157);
}
void car_h_28(double *state, double *unused, double *out_5938536161097710960) {
  h_28(state, unused, out_5938536161097710960);
}
void car_H_28(double *state, double *unused, double *out_1495157348996847574) {
  H_28(state, unused, out_1495157348996847574);
}
void car_h_31(double *state, double *unused, double *out_1792015129981750296) {
  h_31(state, unused, out_1792015129981750296);
}
void car_H_31(double *state, double *unused, double *out_2169904904875118843) {
  H_31(state, unused, out_2169904904875118843);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_init(car);
