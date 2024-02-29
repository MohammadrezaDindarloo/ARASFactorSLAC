#include "scampi_fg_data_types.h"

#pragma once

#include "scampi_fg_IK.h"
#include "scampi_fg_FK.h"
#include "../include/scampi_fg_utils.h"


#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <manif/manif.h>
#include <ceres/ceres.h>

using namespace Eigen;
using namespace std;


void fkSolver(double *lc_cat, 
              Vector3d pos_init,  
              Matrix3d rot_init, Vector2d fc0, 
              RobotParameters<double> params,
              FKDataOut<double> *result);

void ikSolver(RobotParameters<double> params, 
              Eigen::Matrix3d rot_init, 
              Eigen::Vector3d p_platform,
              IKDataOut<double> *result);

void forward_kinematic_factor_graph_optimizer(double lc0, double lc1, double lc2, double lc3,
                              double rot_init_x, double rot_init_y, double rot_init_z, double rot_init_w, 
                              double init_h1, double init_v1, double init_rx,  double init_ry,  double init_rz,
                              double init_tx, double init_ty, double init_tz, gtsam::Values *result_LM);

void inverse_kinematic_factor_graph_optimizer(double p_init_0, double p_init_1, double p_init_2,
                              double rot_init_x, double rot_init_y, double rot_init_z, double rot_init_w, 
                              int largest_cable,
                              double init_h1, double init_v1, double init_rx,  double init_ry,  double init_rz, gtsam::Values *result_LM);
