#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_180823419254858442);
void live_err_fun(double *nom_x, double *delta_x, double *out_382937832037016470);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_6356608129355425497);
void live_H_mod_fun(double *state, double *out_6159620833474390477);
void live_f_fun(double *state, double dt, double *out_6910680473094794463);
void live_F_fun(double *state, double dt, double *out_2569738049418641073);
void live_h_4(double *state, double *unused, double *out_6365147405598968614);
void live_H_4(double *state, double *unused, double *out_8422102642620680189);
void live_h_9(double *state, double *unused, double *out_4866059614826082472);
void live_H_9(double *state, double *unused, double *out_8663292289250270834);
void live_h_10(double *state, double *unused, double *out_5112365462136489217);
void live_H_10(double *state, double *unused, double *out_5631211871600410720);
void live_h_12(double *state, double *unused, double *out_6804380459216229846);
void live_H_12(double *state, double *unused, double *out_5005185023056909632);
void live_h_31(double *state, double *unused, double *out_1382977200322234885);
void live_H_31(double *state, double *unused, double *out_2259621990731895923);
void live_h_32(double *state, double *unused, double *out_8224451516702920237);
void live_H_32(double *state, double *unused, double *out_7570065462321378495);
void live_h_13(double *state, double *unused, double *out_4454733693100208891);
void live_H_13(double *state, double *unused, double *out_5638774433905861398);
void live_h_14(double *state, double *unused, double *out_4866059614826082472);
void live_H_14(double *state, double *unused, double *out_8663292289250270834);
void live_h_33(double *state, double *unused, double *out_8845191609115125344);
void live_H_33(double *state, double *unused, double *out_3507422369077406447);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}