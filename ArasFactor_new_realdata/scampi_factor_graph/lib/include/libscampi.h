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

void forward_kinematic_factor_graph_optimizer(std::vector<Eigen::Matrix<double, 4, 1>> cable_length_collection,
                                            std::vector<Eigen::Matrix<double, 4, 1>> IK_cable_length_collection,
                                            std::vector<Eigen::Matrix<double, 2, 1>> cable_forces_collection,
                                            std::vector<Eigen::Matrix<double, 3, 1>> p_platform_collection,
                                            std::vector<Eigen::Matrix<double, 3, 3>> rot_init_platform_collection,
                                            std::vector<Eigen::Matrix<double, 3, 3>> delta_rot_platform_collection,
                                            Eigen::Matrix<double, 4, 3> pulley_position_estimate,
                                            std::vector<std::vector< Eigen::Matrix<double, 3, 1> >> b_in_w_collection,
                                            std::vector<std::vector< Eigen::Matrix<double, 3, 1> >> p_in_w_collection,
                                            double *optimizer_error,
                                            std::vector<gtsam::Pose3> *Optimized_pose,
                                            std::vector<gtsam::Pose3> *GT_pose,
                                            gtsam::Values *oprimization_result_LM);

void inverse_kinematic_factor_graph_optimizer(Eigen::Vector3d p_init, Eigen::Matrix3d rot_init, int largest_cable,
                              double init_estimate_h1, double init_estimate_v1, gtsam::Rot3 init_estimate_rot, gtsam::Values *oprimization_result_LM);
