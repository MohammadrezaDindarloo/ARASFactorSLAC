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
                              Eigen::Matrix3d rot_init, 
                              double init_estimate_h1, double init_estimate_v1, gtsam::Rot3 init_estimate_rotation,
                              gtsam::Point3 init_estimate_position, gtsam::Values *oprimization_result_LM);

void inverse_kinematic_factor_graph_optimizer(Eigen::Vector3d p_init, Eigen::Matrix3d rot_init, int largest_cable,
                              double init_estimate_h1, double init_estimate_v1, gtsam::Rot3 init_estimate_rot, gtsam::Values *oprimization_result_LM);
