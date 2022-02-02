#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_5096038058135125524);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_8373756370158955549);
void car_H_mod_fun(double *state, double *out_7930719249271816004);
void car_f_fun(double *state, double dt, double *out_2019501809289876349);
void car_F_fun(double *state, double dt, double *out_3574153591255788280);
void car_h_25(double *state, double *unused, double *out_1634828461630621151);
void car_H_25(double *state, double *unused, double *out_2197806516232288857);
void car_h_24(double *state, double *unused, double *out_8272163805170178838);
void car_H_24(double *state, double *unused, double *out_25156917226789291);
void car_h_30(double *state, double *unused, double *out_6049964250792746476);
void car_H_30(double *state, double *unused, double *out_4716139474739537484);
void car_h_26(double *state, double *unused, double *out_4192496911818067350);
void car_H_26(double *state, double *unused, double *out_1543696802641767367);
void car_h_27(double *state, double *unused, double *out_8015606147812933325);
void car_H_27(double *state, double *unused, double *out_4504653125695744252);
void car_h_29(double *state, double *unused, double *out_8172792816164062470);
void car_H_29(double *state, double *unused, double *out_1819658469580927157);
void car_h_28(double *state, double *unused, double *out_5938536161097710960);
void car_H_28(double *state, double *unused, double *out_1495157348996847574);
void car_h_31(double *state, double *unused, double *out_1792015129981750296);
void car_H_31(double *state, double *unused, double *out_2169904904875118843);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}