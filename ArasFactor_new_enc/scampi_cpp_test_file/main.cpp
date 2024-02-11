// Copyright (C) 2021 Rooholla Khorram Bakht, Eren Allak,
// and others, Control of Networked Systems, University of Klagenfurt, Austria.
//
// All rights reserved.
//
// This software is licensed under the terms of the BSD-2-Clause-License with
// no commercial use allowed, the full terms of which are made available
// in the LICENSE file. No license in patents is granted.
//
// You can contact the authors at <r.khorrambakht@alumni.kntu.ac.ir>,
// <eren.allak@aau.at>

// g++ -std=c++17 main.cpp libscampi_ks_solvers.cpp -I/usr/local/include/eigen3 -I/usr/local/include/ceres -L/usr/local/include/ceres -lceres -lcholmod -lstdc++ -lm -lglog -pthread -lmetis -llapack && ./a.out 

#include <ceres/ceres.h>
#include <Eigen/Dense>
#include "libscampi_ks_solvers.h"


using namespace ceres;
using namespace Eigen;

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

class CableRobotParams
{

  public:
    CableRobotParams(const double g_c, const double f_g): f_g_(f_g), g_c_(g_c) {}

    Vector3d p1_, p2_, p3_, p4_;
    Vector3d b1_, b2_, b3_, b4_;
    Vector3d r_to_cog_;
    
    double g_c_, f_g_;
    void setPulleyPoses(Vector3d p1, Vector3d p2, Vector3d p3, Vector3d p4)
    {
      p1_ = p1;
      p2_ = p2;
      p3_ = p3;
      p4_ = p4;
    }
    void setEEAnchors(Vector3d b1, Vector3d b2, Vector3d b3, Vector3d b4)
    {
      b1_ = b1;
      b2_ = b2;
      b3_ = b3;
      b4_ = b4;
    }
    void setCog(Vector3d r_to_cog)
    {
      r_to_cog_ = r_to_cog;
    }
};

std::vector<MatrixXd> inverseKinematicsSolver(CableRobotParams robot_params, Vector3d p_platform, Matrix3d rot_init)
{
  std::vector<MatrixXd> results_list;
  RobotParameters<double> params_;
  //initilize the pulley locations
  params_.pulleys.push_back(robot_params.p1_);
  params_.pulleys.push_back(robot_params.p2_);
  params_.pulleys.push_back(robot_params.p3_);
  params_.pulleys.push_back(robot_params.p4_);

  //initilize cable attachement points
  params_.ef_points.push_back(robot_params.b1_);
  params_.ef_points.push_back(robot_params.b2_);
  params_.ef_points.push_back(robot_params.b3_);
  params_.ef_points.push_back(robot_params.b4_);
  params_.r_to_cog = robot_params.r_to_cog_;
  params_.g_c= robot_params.g_c_;
  params_.f_g = robot_params.f_g_;

  IKDataOut<double> ik_result;
  //Run the first solver to get the end-effector orientation, cable forces, and cable lengths
  ikSolver(p_platform, rot_init, params_, &ik_result);
  //Extract The results and return them as a list of matrices to python
  results_list.push_back(ik_result.rot_platform);
  //Catenary cable lengths
  
  int N = 4;

  VectorXd lc_cat(N);
  lc_cat << ik_result.lc_cat[0], ik_result.lc_cat[1], ik_result.lc_cat[2], ik_result.lc_cat[3];
  results_list.push_back(lc_cat);

  //Cable Forces
  MatrixXd cable_forces(2,N);
  for(int i=0; i<N; i++)
    cable_forces.col(i) =  ik_result.cable_forces[i];
  results_list.push_back(cable_forces);

  //Catenary parameters
  //C1
  MatrixXd c1_(1, N);
  for(int i=0; i<N; i++)
    c1_(0,i) = ik_result.c1[i];
  results_list.push_back(c1_);

  //C2
  MatrixXd c2_(1, N);
  for(int i=0; i<N; i++)
    c2_(0,i) = ik_result.c2[i];
  results_list.push_back(c2_);

  //body coordinates represented in the world frame
  MatrixXd b_in_w_(3,N);
  for(int i=0; i<N; i++)
    b_in_w_.col(i) =  ik_result.b_in_w[i];
  results_list.push_back(b_in_w_);

  return results_list;
}

std::vector<MatrixXd> forwardKinematicsSolver(CableRobotParams robot_params, VectorXd lc_cat, Vector2d fc_1, Vector3d pos_init, Matrix3d rot_init)
{
  std::vector<MatrixXd> results_list;
  RobotParameters<double> params_;
  //initilize the pulley locations
  params_.pulleys.push_back(robot_params.p1_);
  params_.pulleys.push_back(robot_params.p2_);
  params_.pulleys.push_back(robot_params.p3_);
  params_.pulleys.push_back(robot_params.p4_);
  params_.r_to_cog = robot_params.r_to_cog_;

  //initilize cable attachement points
  params_.ef_points.push_back(robot_params.b1_);
  params_.ef_points.push_back(robot_params.b2_);
  params_.ef_points.push_back(robot_params.b3_);
  params_.ef_points.push_back(robot_params.b4_);
  params_.g_c= robot_params.g_c_;
  params_.f_g = robot_params.f_g_;
  //Run the solver
  FKDataOut<double> fk_results;
  fkSolver(lc_cat.data(), pos_init, rot_init, fc_1, params_, &fk_results);

  results_list.push_back(fk_results.rot_platform);
  results_list.push_back(fk_results.p_platform);

  int N = 4;

  //Cable Forces
  MatrixXd cable_forces(2,N);
  for(int i=0; i<N; i++)
    cable_forces.col(i) =  fk_results.cable_forces[i];
  results_list.push_back(cable_forces);

  //Catenary parameters
  //C1
  MatrixXd c1_(1, N);
  for(int i=0; i<N; i++)
    c1_(0,i) = fk_results.c1[i];
  results_list.push_back(c1_);

  //C2
  MatrixXd c2_(1, N);
  for(int i=0; i<N; i++)
    c2_(0,i) = fk_results.c2[i];
  results_list.push_back(c2_);

  //body coordinates represented in the world frame
  MatrixXd b_in_w(3,N);
  for(int i=0; i<N; i++)
    b_in_w.col(i) =  fk_results.b_in_w[i];
  results_list.push_back(b_in_w);

  return results_list;
}

int main(int argc, char const *argv[])
{
  CableRobotParams robo_param(0.1034955, 43.164);
  Eigen::Vector3d Pulley_a(-1.9874742031097412, -8.319656372070312, 8.471846580505371);
  Eigen::Vector3d Pulley_b(2.5193355532036756, -8.388501748709967, 8.469020753679201);
  Eigen::Vector3d Pulley_c(2.717151941069913, 4.774436992746004, 8.364108863330584);
  Eigen::Vector3d Pulley_d(-1.7965602546229, 4.832889384134232, 8.370128714520508);
  robo_param.setPulleyPoses(Pulley_a, Pulley_b, Pulley_c, Pulley_d);

  Eigen::Vector3d Ee_a(-0.21 , -0.21 , -0.011);  
  Eigen::Vector3d Ee_b(0.21  , -0.21 , -0.011);
  Eigen::Vector3d Ee_c(0.21  ,  0.21 , -0.011);
  Eigen::Vector3d Ee_d(-0.21 ,  0.21 , -0.011);
  robo_param.setEEAnchors(Ee_a, Ee_b, Ee_c, Ee_d);

  Eigen::Vector3d r_to_cog(0, 0, -0.12);
  robo_param.setCog(r_to_cog);

  Eigen::Vector3d p_platform(3.09173747e-01, -1.83715841e+00,  -50.18367984e+00);
  Eigen::Matrix3d rot_init;
  rot_init << 0.99268615,  0.11337417, -0.04147891,
             -0.11309773,  0.99354347,  0.00895918,
              0.04222684, -0.00420248,  0.99909921; 
  std::vector<MatrixXd> inverse_results = inverseKinematicsSolver(robo_param, p_platform, rot_init);
  std::cout << std::endl << "---------inverse result--------" << std::endl;
  std::cout << std::endl << "rot_platform: " << std::endl << inverse_results[0] << std::endl;
  std::cout << std::endl << "l_cat: " << std::endl << inverse_results[1] << std::endl;
  std::cout << std::endl << "cable_forces: " << std::endl << inverse_results[2] << std::endl;
  std::cout << std::endl << "c1: " << std::endl << inverse_results[3] << std::endl;
  std::cout << std::endl << "c2: " << std::endl << inverse_results[4] << std::endl;
  std::cout << std::endl << "b_in_w: " << std::endl << inverse_results[5] << std::endl;

  Eigen::VectorXd lc_cat = inverse_results[1];
  Eigen::Vector2d fc_1 = inverse_results[2].col(0);
  Eigen::Vector3d pos_init = p_platform;
  Eigen::Matrix3d rtation_init = rot_init;
  std::vector<MatrixXd> forward_result = forwardKinematicsSolver(robo_param, lc_cat, fc_1, pos_init, rtation_init);

  // std::vector<MatrixXd> forward_result = forwardKinematicsSolver(robo_param, lc_cat, fc_1, pos_init, rtation_init);
  std::cout << std::endl << "---------froward result--------"  << std::endl;
  std::cout << std::endl << "rot_platform: " << std::endl << forward_result[0] << std::endl;
  std::cout << std::endl << "p_platform: " << std::endl << forward_result[1] << std::endl;
  std::cout << std::endl << "cable_forces: " << std::endl << forward_result[2] << std::endl;
  std::cout << std::endl << "c1: " << std::endl << forward_result[3] << std::endl;
  std::cout << std::endl << "c2: " << std::endl << forward_result[4] << std::endl;
  std::cout << std::endl << "b_in_w: " << std::endl << forward_result[5] << std::endl;
  std::cout << std::endl << "p_platform_error in  mm: " << std::endl << (forward_result[1] - p_platform).norm() * 1000 << std::endl;

  return 0;
}
