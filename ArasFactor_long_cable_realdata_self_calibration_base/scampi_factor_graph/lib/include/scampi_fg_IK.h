#pragma once

#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <boost/optional.hpp>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/inference/Symbol.h>
#include "scampi_function_header_include.h"
#include <gtsam/base/Matrix.h>
#include <gtsam/geometry/Pose3.h>
#include <gtsam/geometry/Rot3.h>
#include <gtsam/slam/PriorFactor.h>
#include <gtsam/slam/BetweenFactor.h>
#include <sym/rot3.h>
#include <sym/pose3.h>

using namespace gtsam;
using namespace std;

Eigen::Matrix3d gtsamRot3ToEigenMatrix(const gtsam::Rot3& rot3) {
    Eigen::Matrix3d eigenMatrix;
    eigenMatrix << rot3.matrix();
    return eigenMatrix;
}

sym::Rot3<double> SymforceFromGtsam(const gtsam::Rot3& gtsam_rot3) {
  return sym::Rot3<double>(gtsam_rot3.toQuaternion());
}

// Function to convert a GTSAM Pose3 to a symbolic Pose3
sym::Pose3<double> SymPoseFromGtsam(const gtsam::Pose3& gtsam_pose) {
    // Extract rotation matrix and translation vector from GTSAM Pose3
    gtsam::Rot3 rotation = gtsam_pose.rotation();
    gtsam::Point3 translation = gtsam_pose.translation();

    sym::Pose3<double> symPose;

    symPose = sym::Pose3(SymforceFromGtsam(rotation), translation);

    return symPose;
}

gtsam::Rot3 EigenMatrixToGtsamRot3(const Eigen::Matrix<double, 3, 3>& eigenMatrix) {
    // Create a GTSAM Rot3 object from the Eigen matrix
    gtsam::Rot3 gtsamRotation(eigenMatrix);
    
    return gtsamRotation;
}

Eigen::Matrix<double, 4, 3> EigenFromGtsam43(const gtsam::Matrix& gtsam_matrix) {
    // Assuming gtsam::Matrix is a class provided by GTSAM for matrices
    // You may need to adjust this based on the actual class in GTSAM

    // Ensure the GTSAM matrix is 4x3
    if (gtsam_matrix.rows() != 4 || gtsam_matrix.cols() != 3) {
        throw std::invalid_argument("Input matrix must be 4x3");
    }

    // Map the GTSAM matrix to Eigen matrix
    Eigen::Map<const Eigen::Matrix<double, 4, 3>> eigen_matrix(gtsam_matrix.data(), 4, 3);

    return eigen_matrix;
}


namespace gtsam
{
    class IK_factor_graoh_cost1 : public NoiseModelFactorN<double, double, gtsam::Rot3, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Point3>
    {

    private:
        gtsam::Point3 p_init;
        Eigen::Matrix3d rot_init;
        int largest_cable = 0;

    public:
        // Constructor
        IK_factor_graoh_cost1(Key key1, Key key2, Key key3, Key key4, Key key5, Key key6, Key key7, gtsam::Point3 p_init_, Eigen::Matrix3d rot_init_, const int largest_cable_, const SharedNoiseModel &model) 
        : NoiseModelFactorN<double, double, gtsam::Rot3, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Point3>(model, key1, key2, key3, key4, key5, key6, key7), p_init(p_init_), rot_init(rot_init_), largest_cable(largest_cable_) {}
        
        // Evaluate the error
        Vector evaluateError(const double &fh1, const double &fv1, const gtsam::Rot3 &DeltaRot, const gtsam::Point3 &p_a, const gtsam::Point3 &p_b, const gtsam::Point3 &p_c, const gtsam::Point3 &p_d, 
                             OptionalMatrixType H1,
                             OptionalMatrixType H2,
                             OptionalMatrixType H3,
                             OptionalMatrixType H4,
                             OptionalMatrixType H5,
                             OptionalMatrixType H6,
                             OptionalMatrixType H7) const override
        {   
            Eigen::Matrix<double, 4, 1> Ikresidual_func;

            switch (largest_cable)

            {   
                case 0:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl0 = sym::IkResidualFuncCost1WrtFh1Nl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl0).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl0 = sym::IkResidualFuncCost1WrtFv1Nl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl0).finished();
                    }
                    if (H3) 
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl0 = sym::IkResidualFuncCost1WrtDeltarotNl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl0).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl0 = sym::IkResidualFuncCost1WrtPaNl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl0).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl0 = sym::IkResidualFuncCost1WrtPbNl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl0).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl0 = sym::IkResidualFuncCost1WrtPcNl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl0).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl0 = sym::IkResidualFuncCost1WrtPdNl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl0).finished();
                    }
                    break;


                case 1:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl1 = sym::IkResidualFuncCost1WrtFh1Nl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl1).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl1 = sym::IkResidualFuncCost1WrtFv1Nl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl1).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl1 = sym::IkResidualFuncCost1WrtDeltarotNl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl1).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl1 = sym::IkResidualFuncCost1WrtPaNl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl1).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl1 = sym::IkResidualFuncCost1WrtPbNl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl1).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl1 = sym::IkResidualFuncCost1WrtPcNl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl1).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl1 = sym::IkResidualFuncCost1WrtPdNl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl1).finished();
                    }
                    break;


                case 2:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl2 = sym::IkResidualFuncCost1WrtFh1Nl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl2).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl2 = sym::IkResidualFuncCost1WrtFv1Nl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl2).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl2 = sym::IkResidualFuncCost1WrtDeltarotNl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl2).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl2 = sym::IkResidualFuncCost1WrtPaNl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl2).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl2 = sym::IkResidualFuncCost1WrtPbNl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl2).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl2 = sym::IkResidualFuncCost1WrtPcNl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl2).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl2 = sym::IkResidualFuncCost1WrtPdNl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl2).finished();
                    }
                    break;



                case 3:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl3(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl3 = sym::IkResidualFuncCost1WrtFh1Nl3(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl3).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl3 = sym::IkResidualFuncCost1WrtFv1Nl3(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl3).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl3 = sym::IkResidualFuncCost1WrtDeltarotNl3(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl3).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl3 = sym::IkResidualFuncCost1WrtPaNl3(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl3).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl3 = sym::IkResidualFuncCost1WrtPbNl3(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl3).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl3 = sym::IkResidualFuncCost1WrtPcNl3(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl3).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl3 = sym::IkResidualFuncCost1WrtPdNl3(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl3).finished();
                    }
                    break;



                case 4:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl4 = sym::IkResidualFuncCost1WrtFh1Nl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl4).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl4 = sym::IkResidualFuncCost1WrtFv1Nl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl4).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl4 = sym::IkResidualFuncCost1WrtDeltarotNl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl4).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl4 = sym::IkResidualFuncCost1WrtPaNl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl4).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl4 = sym::IkResidualFuncCost1WrtPbNl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl4).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl4 = sym::IkResidualFuncCost1WrtPcNl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl4).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl4 = sym::IkResidualFuncCost1WrtPdNl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl4).finished();
                    }
                    break;


                case 5:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl5(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl5 = sym::IkResidualFuncCost1WrtFh1Nl5(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl5).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl5 = sym::IkResidualFuncCost1WrtFv1Nl5(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl5).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl5 = sym::IkResidualFuncCost1WrtDeltarotNl5(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl5).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl5 = sym::IkResidualFuncCost1WrtPaNl5(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl5).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl5 = sym::IkResidualFuncCost1WrtPbNl5(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl5).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl5 = sym::IkResidualFuncCost1WrtPcNl5(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl5).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl5 = sym::IkResidualFuncCost1WrtPdNl5(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl5).finished();
                    }
                    break;


                case 6:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl6 = sym::IkResidualFuncCost1WrtFh1Nl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl6).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl6 = sym::IkResidualFuncCost1WrtFv1Nl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl6).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl6 = sym::IkResidualFuncCost1WrtDeltarotNl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl6).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl6 = sym::IkResidualFuncCost1WrtPaNl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl6).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl6 = sym::IkResidualFuncCost1WrtPbNl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl6).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl6 = sym::IkResidualFuncCost1WrtPcNl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl6).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl6 = sym::IkResidualFuncCost1WrtPdNl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl6).finished();
                    }
                    break;



                case 7:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl7(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl7 = sym::IkResidualFuncCost1WrtFh1Nl7(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl7).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl7 = sym::IkResidualFuncCost1WrtFv1Nl7(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl7).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl7 = sym::IkResidualFuncCost1WrtDeltarotNl7(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl7).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl7 = sym::IkResidualFuncCost1WrtPaNl7(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl7).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl7 = sym::IkResidualFuncCost1WrtPbNl7(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl7).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl7 = sym::IkResidualFuncCost1WrtPcNl7(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl7).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl7 = sym::IkResidualFuncCost1WrtPdNl7(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl7).finished();
                    }
                    break;


                case 8:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl8 = sym::IkResidualFuncCost1WrtFh1Nl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl8).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl8 = sym::IkResidualFuncCost1WrtFv1Nl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl8).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl8 = sym::IkResidualFuncCost1WrtDeltarotNl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl8).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl8 = sym::IkResidualFuncCost1WrtPaNl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl8).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl8 = sym::IkResidualFuncCost1WrtPbNl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl8).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl8 = sym::IkResidualFuncCost1WrtPcNl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl8).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl8 = sym::IkResidualFuncCost1WrtPdNl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl8).finished();
                    }
                    break;



                case 9:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl9 = sym::IkResidualFuncCost1WrtFh1Nl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl9).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl9 = sym::IkResidualFuncCost1WrtFv1Nl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl9).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl9 = sym::IkResidualFuncCost1WrtDeltarotNl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl9).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl9 = sym::IkResidualFuncCost1WrtPaNl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl9).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl9 = sym::IkResidualFuncCost1WrtPbNl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl9).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl9 = sym::IkResidualFuncCost1WrtPcNl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl9).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl9 = sym::IkResidualFuncCost1WrtPdNl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl9).finished();
                    }
                    break;



                case 10:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl10 = sym::IkResidualFuncCost1WrtFh1Nl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl10).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl10 = sym::IkResidualFuncCost1WrtFv1Nl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl10).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl10 = sym::IkResidualFuncCost1WrtDeltarotNl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl10).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl10 = sym::IkResidualFuncCost1WrtPaNl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl10).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl10 = sym::IkResidualFuncCost1WrtPbNl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl10).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl10 = sym::IkResidualFuncCost1WrtPcNl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl10).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl10 = sym::IkResidualFuncCost1WrtPdNl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl10).finished();
                    }
                    break;



                case 11:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl11 = sym::IkResidualFuncCost1WrtFh1Nl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl11).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl11 = sym::IkResidualFuncCost1WrtFv1Nl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl11).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl11 = sym::IkResidualFuncCost1WrtDeltarotNl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl11).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl11 = sym::IkResidualFuncCost1WrtPaNl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl11).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl11 = sym::IkResidualFuncCost1WrtPbNl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl11).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl11 = sym::IkResidualFuncCost1WrtPcNl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl11).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl11 = sym::IkResidualFuncCost1WrtPdNl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl11).finished();
                    }
                    break;



                case 12:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl12 = sym::IkResidualFuncCost1WrtFh1Nl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl12).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl12 = sym::IkResidualFuncCost1WrtFv1Nl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl12).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl12 = sym::IkResidualFuncCost1WrtDeltarotNl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl12).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl12 = sym::IkResidualFuncCost1WrtPaNl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl12).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl12 = sym::IkResidualFuncCost1WrtPbNl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl12).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl12 = sym::IkResidualFuncCost1WrtPcNl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl12).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl12 = sym::IkResidualFuncCost1WrtPdNl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl12).finished();
                    }
                    break;



                case 13:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl13 = sym::IkResidualFuncCost1WrtFh1Nl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl13).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl13 = sym::IkResidualFuncCost1WrtFv1Nl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl13).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl13 = sym::IkResidualFuncCost1WrtDeltarotNl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl13).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl13 = sym::IkResidualFuncCost1WrtPaNl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl13).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl13 = sym::IkResidualFuncCost1WrtPbNl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl13).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl13 = sym::IkResidualFuncCost1WrtPcNl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl13).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl13 = sym::IkResidualFuncCost1WrtPdNl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl13).finished();
                    }
                    break;



                case 14:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl14 = sym::IkResidualFuncCost1WrtFh1Nl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl14).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl14 = sym::IkResidualFuncCost1WrtFv1Nl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl14).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl14 = sym::IkResidualFuncCost1WrtDeltarotNl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl14).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl14 = sym::IkResidualFuncCost1WrtPaNl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl14).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl14 = sym::IkResidualFuncCost1WrtPbNl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl14).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl14 = sym::IkResidualFuncCost1WrtPcNl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl14).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl14 = sym::IkResidualFuncCost1WrtPdNl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl14).finished();
                    }
                    break;



                case 15:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl15 = sym::IkResidualFuncCost1WrtFh1Nl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl15).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl15 = sym::IkResidualFuncCost1WrtFv1Nl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl15).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl15 = sym::IkResidualFuncCost1WrtDeltarotNl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl15).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl15 = sym::IkResidualFuncCost1WrtPaNl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl15).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl15 = sym::IkResidualFuncCost1WrtPbNl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl15).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl15 = sym::IkResidualFuncCost1WrtPcNl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl15).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl15 = sym::IkResidualFuncCost1WrtPdNl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl15).finished();
                    }
                    break;



                case 16:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl16 = sym::IkResidualFuncCost1WrtFh1Nl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl16).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl16 = sym::IkResidualFuncCost1WrtFv1Nl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl16).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl16 = sym::IkResidualFuncCost1WrtDeltarotNl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl16).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl16 = sym::IkResidualFuncCost1WrtPaNl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl16).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl16 = sym::IkResidualFuncCost1WrtPbNl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl16).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl16 = sym::IkResidualFuncCost1WrtPcNl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl16).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl16 = sym::IkResidualFuncCost1WrtPdNl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl16).finished();
                    }
                    break;



                case 17:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl17 = sym::IkResidualFuncCost1WrtFh1Nl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl17).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl17 = sym::IkResidualFuncCost1WrtFv1Nl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl17).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl17 = sym::IkResidualFuncCost1WrtDeltarotNl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl17).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl17 = sym::IkResidualFuncCost1WrtPaNl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl17).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl17 = sym::IkResidualFuncCost1WrtPbNl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl17).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl17 = sym::IkResidualFuncCost1WrtPcNl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl17).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl17 = sym::IkResidualFuncCost1WrtPdNl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl17).finished();
                    }
                    break;


                case 18:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl18 = sym::IkResidualFuncCost1WrtFh1Nl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl18).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl18 = sym::IkResidualFuncCost1WrtFv1Nl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl18).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl18 = sym::IkResidualFuncCost1WrtDeltarotNl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl18).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl18 = sym::IkResidualFuncCost1WrtPaNl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl18).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl18 = sym::IkResidualFuncCost1WrtPbNl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl18).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl18 = sym::IkResidualFuncCost1WrtPcNl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl18).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl18 = sym::IkResidualFuncCost1WrtPdNl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl18).finished();
                    }
                    break;



                case 19:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl19 = sym::IkResidualFuncCost1WrtFh1Nl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl19).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl19 = sym::IkResidualFuncCost1WrtFv1Nl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl19).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl19 = sym::IkResidualFuncCost1WrtDeltarotNl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl19).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl19 = sym::IkResidualFuncCost1WrtPaNl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl19).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl19 = sym::IkResidualFuncCost1WrtPbNl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl19).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl19 = sym::IkResidualFuncCost1WrtPcNl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl19).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl19 = sym::IkResidualFuncCost1WrtPdNl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl19).finished();
                    }
                    break;



                case 20:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl20 = sym::IkResidualFuncCost1WrtFh1Nl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl20).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl20 = sym::IkResidualFuncCost1WrtFv1Nl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl20).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl20 = sym::IkResidualFuncCost1WrtDeltarotNl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl20).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl20 = sym::IkResidualFuncCost1WrtPaNl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl20).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl20 = sym::IkResidualFuncCost1WrtPbNl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl20).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl20 = sym::IkResidualFuncCost1WrtPcNl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl20).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl20 = sym::IkResidualFuncCost1WrtPdNl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl20).finished();
                    }
                    break;



                case 21:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl21 = sym::IkResidualFuncCost1WrtFh1Nl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl21).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl21 = sym::IkResidualFuncCost1WrtFv1Nl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl21).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl21 = sym::IkResidualFuncCost1WrtDeltarotNl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl21).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl21 = sym::IkResidualFuncCost1WrtPaNl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl21).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl21 = sym::IkResidualFuncCost1WrtPbNl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl21).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl21 = sym::IkResidualFuncCost1WrtPcNl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl21).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl21 = sym::IkResidualFuncCost1WrtPdNl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl21).finished();
                    }
                    break;



                case 22:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl22 = sym::IkResidualFuncCost1WrtFh1Nl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl22).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl22 = sym::IkResidualFuncCost1WrtFv1Nl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl22).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl22 = sym::IkResidualFuncCost1WrtDeltarotNl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl22).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl22 = sym::IkResidualFuncCost1WrtPaNl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl22).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl22 = sym::IkResidualFuncCost1WrtPbNl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl22).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl22 = sym::IkResidualFuncCost1WrtPcNl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl22).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl22 = sym::IkResidualFuncCost1WrtPdNl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl22).finished();
                    }
                    break;



                case 23:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl23 = sym::IkResidualFuncCost1WrtFh1Nl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl23).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl23 = sym::IkResidualFuncCost1WrtFv1Nl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl23).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl23 = sym::IkResidualFuncCost1WrtDeltarotNl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl23).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl23 = sym::IkResidualFuncCost1WrtPaNl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl23).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl23 = sym::IkResidualFuncCost1WrtPbNl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl23).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl23 = sym::IkResidualFuncCost1WrtPcNl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl23).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl23 = sym::IkResidualFuncCost1WrtPdNl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl23).finished();
                    }
                    break;
            }
            // std::cout << "cost1: " << Ikresidual_func << std::endl;
            return (Vector(4) << Ikresidual_func).finished();
        }
    };

    class IK_factor_graoh_cost2 : public NoiseModelFactorN<double, double, gtsam::Rot3, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Point3>
    {

    private:
        gtsam::Point3 p_init;
        Eigen::Matrix3d rot_init;
        int largest_cable = 0;

    public:
        // Constructor
        IK_factor_graoh_cost2(Key key1, Key key2, Key key3, Key key4, Key key5, Key key6, Key key7, gtsam::Point3 p_init_, Eigen::Matrix3d rot_init_, const int largest_cable_, const SharedNoiseModel &model) 
        : NoiseModelFactorN<double, double, gtsam::Rot3, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Point3>(model, key1, key2, key3, key4, key5, key6, key7), p_init(p_init_), rot_init(rot_init_), largest_cable(largest_cable_) {}

        // Evaluate the error
        Vector evaluateError(const double &fh1, const double &fv1, const gtsam::Rot3 &DeltaRot, const gtsam::Point3 &p_a, const gtsam::Point3 &p_b, const gtsam::Point3 &p_c, const gtsam::Point3 &p_d,
                             OptionalMatrixType H1,
                             OptionalMatrixType H2,
                             OptionalMatrixType H3,
                             OptionalMatrixType H4,
                             OptionalMatrixType H5,
                             OptionalMatrixType H6,
                             OptionalMatrixType H7) const override
        {   
            Eigen::Matrix<double, 4, 1> Ikresidual_func;

            switch (largest_cable)
            {   
                case 0:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl0 = sym::IkResidualFuncCost2WrtFh1Nl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl0).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl0 = sym::IkResidualFuncCost2WrtFv1Nl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl0).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl0 = sym::IkResidualFuncCost2WrtDeltarotNl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl0).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl0 = sym::IkResidualFuncCost2WrtPaNl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl0).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl0 = sym::IkResidualFuncCost2WrtPbNl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl0).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl0 = sym::IkResidualFuncCost2WrtPcNl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl0).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl0 = sym::IkResidualFuncCost2WrtPdNl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl0).finished();
                    }
                    break;


                case 1:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl1 = sym::IkResidualFuncCost2WrtFh1Nl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl1).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl1 = sym::IkResidualFuncCost2WrtFv1Nl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl1).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl1 = sym::IkResidualFuncCost2WrtDeltarotNl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl1).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl1 = sym::IkResidualFuncCost2WrtPaNl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl1).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl1 = sym::IkResidualFuncCost2WrtPbNl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl1).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl1 = sym::IkResidualFuncCost2WrtPcNl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl1).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl1 = sym::IkResidualFuncCost2WrtPdNl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl1).finished();
                    }
                    break;


                case 2:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl2 = sym::IkResidualFuncCost2WrtFh1Nl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl2).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl2 = sym::IkResidualFuncCost2WrtFv1Nl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl2).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl2 = sym::IkResidualFuncCost2WrtDeltarotNl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl2).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl2 = sym::IkResidualFuncCost2WrtPaNl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl2).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl2 = sym::IkResidualFuncCost2WrtPbNl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl2).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl2 = sym::IkResidualFuncCost2WrtPcNl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl2).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl2 = sym::IkResidualFuncCost2WrtPdNl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl2).finished();
                    }
                    break;



                case 3:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl3(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl3 = sym::IkResidualFuncCost2WrtFh1Nl3(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl3).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl3 = sym::IkResidualFuncCost2WrtFv1Nl3(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl3).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl3 = sym::IkResidualFuncCost2WrtDeltarotNl3(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl3).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl3 = sym::IkResidualFuncCost2WrtPaNl3(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl3).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl3 = sym::IkResidualFuncCost2WrtPbNl3(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl3).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl3 = sym::IkResidualFuncCost2WrtPcNl3(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl3).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl3 = sym::IkResidualFuncCost2WrtPdNl3(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl3).finished();
                    }
                    break;



                case 4:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl4 = sym::IkResidualFuncCost2WrtFh1Nl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl4).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl4 = sym::IkResidualFuncCost2WrtFv1Nl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl4).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl4 = sym::IkResidualFuncCost2WrtDeltarotNl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl4).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl4 = sym::IkResidualFuncCost2WrtPaNl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl4).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl4 = sym::IkResidualFuncCost2WrtPbNl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl4).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl4 = sym::IkResidualFuncCost2WrtPcNl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl4).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl4 = sym::IkResidualFuncCost2WrtPdNl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl4).finished();
                    }
                    break;


                case 5:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl5(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl5 = sym::IkResidualFuncCost2WrtFh1Nl5(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl5).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl5 = sym::IkResidualFuncCost2WrtFv1Nl5(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl5).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl5 = sym::IkResidualFuncCost2WrtDeltarotNl5(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl5).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl5 = sym::IkResidualFuncCost2WrtPaNl5(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl5).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl5 = sym::IkResidualFuncCost2WrtPbNl5(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl5).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl5 = sym::IkResidualFuncCost2WrtPcNl5(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl5).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl5 = sym::IkResidualFuncCost2WrtPdNl5(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl5).finished();
                    }
                    break;


                case 6:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl6 = sym::IkResidualFuncCost2WrtFh1Nl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl6).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl6 = sym::IkResidualFuncCost2WrtFv1Nl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl6).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl6 = sym::IkResidualFuncCost2WrtDeltarotNl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl6).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl6 = sym::IkResidualFuncCost2WrtPaNl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl6).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl6 = sym::IkResidualFuncCost2WrtPbNl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl6).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl6 = sym::IkResidualFuncCost2WrtPcNl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl6).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl6 = sym::IkResidualFuncCost2WrtPdNl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl6).finished();
                    }
                    break;



                case 7:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl7(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl7 = sym::IkResidualFuncCost2WrtFh1Nl7(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl7).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl7 = sym::IkResidualFuncCost2WrtFv1Nl7(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl7).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl7 = sym::IkResidualFuncCost2WrtDeltarotNl7(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl7).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl7 = sym::IkResidualFuncCost2WrtPaNl7(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl7).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl7 = sym::IkResidualFuncCost2WrtPbNl7(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl7).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl7 = sym::IkResidualFuncCost2WrtPcNl7(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl7).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl7 = sym::IkResidualFuncCost2WrtPdNl7(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl7).finished();
                    }
                    break;


                case 8:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl8 = sym::IkResidualFuncCost2WrtFh1Nl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl8).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl8 = sym::IkResidualFuncCost2WrtFv1Nl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl8).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl8 = sym::IkResidualFuncCost2WrtDeltarotNl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl8).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl8 = sym::IkResidualFuncCost2WrtPaNl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl8).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl8 = sym::IkResidualFuncCost2WrtPbNl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl8).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl8 = sym::IkResidualFuncCost2WrtPcNl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl8).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl8 = sym::IkResidualFuncCost2WrtPdNl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl8).finished();
                    }
                    break;



                case 9:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl9 = sym::IkResidualFuncCost2WrtFh1Nl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl9).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl9 = sym::IkResidualFuncCost2WrtFv1Nl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl9).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl9 = sym::IkResidualFuncCost2WrtDeltarotNl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl9).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl9 = sym::IkResidualFuncCost2WrtPaNl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl9).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl9 = sym::IkResidualFuncCost2WrtPbNl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl9).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl9 = sym::IkResidualFuncCost2WrtPcNl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl9).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl9 = sym::IkResidualFuncCost2WrtPdNl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl9).finished();
                    }
                    break;



                case 10:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl10 = sym::IkResidualFuncCost2WrtFh1Nl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl10).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl10 = sym::IkResidualFuncCost2WrtFv1Nl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl10).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl10 = sym::IkResidualFuncCost2WrtDeltarotNl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl10).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl10 = sym::IkResidualFuncCost2WrtPaNl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl10).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl10 = sym::IkResidualFuncCost2WrtPbNl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl10).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl10 = sym::IkResidualFuncCost2WrtPcNl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl10).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl10 = sym::IkResidualFuncCost2WrtPdNl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl10).finished();
                    }
                    break;



                case 11:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl11 = sym::IkResidualFuncCost2WrtFh1Nl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl11).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl11 = sym::IkResidualFuncCost2WrtFv1Nl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl11).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl11 = sym::IkResidualFuncCost2WrtDeltarotNl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl11).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl11 = sym::IkResidualFuncCost2WrtPaNl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl11).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl11 = sym::IkResidualFuncCost2WrtPbNl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl11).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl11 = sym::IkResidualFuncCost2WrtPcNl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl11).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl11 = sym::IkResidualFuncCost2WrtPdNl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl11).finished();
                    }
                    break;



                case 12:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl12 = sym::IkResidualFuncCost2WrtFh1Nl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl12).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl12 = sym::IkResidualFuncCost2WrtFv1Nl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl12).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl12 = sym::IkResidualFuncCost2WrtDeltarotNl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl12).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl12 = sym::IkResidualFuncCost2WrtPaNl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl12).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl12 = sym::IkResidualFuncCost2WrtPbNl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl12).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl12 = sym::IkResidualFuncCost2WrtPcNl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl12).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl12 = sym::IkResidualFuncCost2WrtPdNl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl12).finished();
                    }
                    break;



                case 13:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl13 = sym::IkResidualFuncCost2WrtFh1Nl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl13).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl13 = sym::IkResidualFuncCost2WrtFv1Nl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl13).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl13 = sym::IkResidualFuncCost2WrtDeltarotNl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl13).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl13 = sym::IkResidualFuncCost2WrtPaNl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl13).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl13 = sym::IkResidualFuncCost2WrtPbNl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl13).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl13 = sym::IkResidualFuncCost2WrtPcNl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl13).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl13 = sym::IkResidualFuncCost2WrtPdNl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl13).finished();
                    }
                    break;



                case 14:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl14 = sym::IkResidualFuncCost2WrtFh1Nl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl14).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl14 = sym::IkResidualFuncCost2WrtFv1Nl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl14).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl14 = sym::IkResidualFuncCost2WrtDeltarotNl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl14).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl14 = sym::IkResidualFuncCost2WrtPaNl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl14).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl14 = sym::IkResidualFuncCost2WrtPbNl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl14).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl14 = sym::IkResidualFuncCost2WrtPcNl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl14).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl14 = sym::IkResidualFuncCost2WrtPdNl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl14).finished();
                    }
                    break;



                case 15:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl15 = sym::IkResidualFuncCost2WrtFh1Nl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl15).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl15 = sym::IkResidualFuncCost2WrtFv1Nl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl15).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl15 = sym::IkResidualFuncCost2WrtDeltarotNl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl15).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl15 = sym::IkResidualFuncCost2WrtPaNl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl15).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl15 = sym::IkResidualFuncCost2WrtPbNl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl15).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl15 = sym::IkResidualFuncCost2WrtPcNl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl15).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl15 = sym::IkResidualFuncCost2WrtPdNl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl15).finished();
                    }
                    break;



                case 16:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl16 = sym::IkResidualFuncCost2WrtFh1Nl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl16).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl16 = sym::IkResidualFuncCost2WrtFv1Nl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl16).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl16 = sym::IkResidualFuncCost2WrtDeltarotNl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl16).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl16 = sym::IkResidualFuncCost2WrtPaNl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl16).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl16 = sym::IkResidualFuncCost2WrtPbNl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl16).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl16 = sym::IkResidualFuncCost2WrtPcNl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl16).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl16 = sym::IkResidualFuncCost2WrtPdNl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl16).finished();
                    }
                    break;



                case 17:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl17 = sym::IkResidualFuncCost2WrtFh1Nl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl17).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl17 = sym::IkResidualFuncCost2WrtFv1Nl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl17).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl17 = sym::IkResidualFuncCost2WrtDeltarotNl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl17).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl17 = sym::IkResidualFuncCost2WrtPaNl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl17).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl17 = sym::IkResidualFuncCost2WrtPbNl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl17).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl17 = sym::IkResidualFuncCost2WrtPcNl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl17).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl17 = sym::IkResidualFuncCost2WrtPdNl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl17).finished();
                    }
                    break;


                case 18:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl18 = sym::IkResidualFuncCost2WrtFh1Nl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl18).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl18 = sym::IkResidualFuncCost2WrtFv1Nl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl18).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl18 = sym::IkResidualFuncCost2WrtDeltarotNl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl18).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl18 = sym::IkResidualFuncCost2WrtPaNl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl18).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl18 = sym::IkResidualFuncCost2WrtPbNl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl18).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl18 = sym::IkResidualFuncCost2WrtPcNl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl18).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl18 = sym::IkResidualFuncCost2WrtPdNl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl18).finished();
                    }
                    break;



                case 19:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl19 = sym::IkResidualFuncCost2WrtFh1Nl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl19).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl19 = sym::IkResidualFuncCost2WrtFv1Nl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl19).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl19 = sym::IkResidualFuncCost2WrtDeltarotNl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl19).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl19 = sym::IkResidualFuncCost2WrtPaNl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl19).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl19 = sym::IkResidualFuncCost2WrtPbNl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl19).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl19 = sym::IkResidualFuncCost2WrtPcNl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl19).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl19 = sym::IkResidualFuncCost2WrtPdNl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl19).finished();
                    }
                    break;



                case 20:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl20 = sym::IkResidualFuncCost2WrtFh1Nl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl20).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl20 = sym::IkResidualFuncCost2WrtFv1Nl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl20).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl20 = sym::IkResidualFuncCost2WrtDeltarotNl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl20).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl20 = sym::IkResidualFuncCost2WrtPaNl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl20).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl20 = sym::IkResidualFuncCost2WrtPbNl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl20).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl20 = sym::IkResidualFuncCost2WrtPcNl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl20).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl20 = sym::IkResidualFuncCost2WrtPdNl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl20).finished();
                    }
                    break;



                case 21:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl21 = sym::IkResidualFuncCost2WrtFh1Nl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl21).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl21 = sym::IkResidualFuncCost2WrtFv1Nl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl21).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl21 = sym::IkResidualFuncCost2WrtDeltarotNl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl21).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl21 = sym::IkResidualFuncCost2WrtPaNl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl21).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl21 = sym::IkResidualFuncCost2WrtPbNl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl21).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl21 = sym::IkResidualFuncCost2WrtPcNl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl21).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl21 = sym::IkResidualFuncCost2WrtPdNl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl21).finished();
                    }
                    break;



                case 22:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl22 = sym::IkResidualFuncCost2WrtFh1Nl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl22).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl22 = sym::IkResidualFuncCost2WrtFv1Nl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl22).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl22 = sym::IkResidualFuncCost2WrtDeltarotNl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl22).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl22 = sym::IkResidualFuncCost2WrtPaNl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl22).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl22 = sym::IkResidualFuncCost2WrtPbNl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl22).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl22 = sym::IkResidualFuncCost2WrtPcNl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl22).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl22 = sym::IkResidualFuncCost2WrtPdNl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl22).finished();
                    }
                    break;



                case 23:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl23 = sym::IkResidualFuncCost2WrtFh1Nl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl23).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl23 = sym::IkResidualFuncCost2WrtFv1Nl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl23).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl23 = sym::IkResidualFuncCost2WrtDeltarotNl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl23).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl23 = sym::IkResidualFuncCost2WrtPaNl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl23).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl23 = sym::IkResidualFuncCost2WrtPbNl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl23).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl23 = sym::IkResidualFuncCost2WrtPcNl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl23).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl23 = sym::IkResidualFuncCost2WrtPdNl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl23).finished();
                    }
                    break;
            }
            // std::cout << "cost2: " << Ikresidual_func << std::endl;
            return (Vector(4) << Ikresidual_func).finished(); 
        }
    };


    class IK_factor_graoh_cost3 : public NoiseModelFactorN<double, double, gtsam::Rot3, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Point3>
    {

    private:
        gtsam::Point3 p_init;
        Eigen::Matrix3d rot_init;
        int largest_cable;

    public:
        // Constructor
        IK_factor_graoh_cost3(Key key1, Key key2, Key key3, Key key4, Key key5, Key key6, Key key7, gtsam::Point3 p_init_, Eigen::Matrix3d rot_init_, const int largest_cable_, const SharedNoiseModel &model) 
        : NoiseModelFactorN<double, double, gtsam::Rot3, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Point3>(model, key1, key2, key3, key4, key5, key6, key7), p_init(p_init_), rot_init(rot_init_), largest_cable(largest_cable_) {}

        // Evaluate the error
        Vector evaluateError(const double &fh1, const double &fv1, const gtsam::Rot3 &DeltaRot, const gtsam::Point3 &p_a, const gtsam::Point3 &p_b, const gtsam::Point3 &p_c, const gtsam::Point3 &p_d,
                             OptionalMatrixType H1,
                             OptionalMatrixType H2,
                             OptionalMatrixType H3,
                             OptionalMatrixType H4,
                             OptionalMatrixType H5,
                             OptionalMatrixType H6,
                             OptionalMatrixType H7) const override
        {   
            Eigen::Matrix<double, 4, 1> Ikresidual_func; 

            switch (largest_cable)
            {   
                case 0:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl0 = sym::IkResidualFuncCost3WrtFh1Nl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl0).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl0 = sym::IkResidualFuncCost3WrtFv1Nl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl0).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl0 = sym::IkResidualFuncCost3WrtDeltarotNl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl0).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl0 = sym::IkResidualFuncCost3WrtPaNl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl0).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl0 = sym::IkResidualFuncCost3WrtPbNl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl0).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl0 = sym::IkResidualFuncCost3WrtPcNl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl0).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl0 = sym::IkResidualFuncCost3WrtPdNl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl0).finished();
                    }
                    break;


                case 1:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl1 = sym::IkResidualFuncCost3WrtFh1Nl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl1).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl1 = sym::IkResidualFuncCost3WrtFv1Nl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl1).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl1 = sym::IkResidualFuncCost3WrtDeltarotNl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl1).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl1 = sym::IkResidualFuncCost3WrtPaNl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl1).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl1 = sym::IkResidualFuncCost3WrtPbNl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl1).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl1 = sym::IkResidualFuncCost3WrtPcNl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl1).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl1 = sym::IkResidualFuncCost3WrtPdNl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl1).finished();
                    }
                    break;


                case 2:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl2 = sym::IkResidualFuncCost3WrtFh1Nl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl2).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl2 = sym::IkResidualFuncCost3WrtFv1Nl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl2).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl2 = sym::IkResidualFuncCost3WrtDeltarotNl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl2).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl2 = sym::IkResidualFuncCost3WrtPaNl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl2).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl2 = sym::IkResidualFuncCost3WrtPbNl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl2).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl2 = sym::IkResidualFuncCost3WrtPcNl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl2).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl2 = sym::IkResidualFuncCost3WrtPdNl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl2).finished();
                    }
                    break;



                case 3:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl3(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl3 = sym::IkResidualFuncCost3WrtFh1Nl3(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl3).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl3 = sym::IkResidualFuncCost3WrtFv1Nl3(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl3).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl3 = sym::IkResidualFuncCost3WrtDeltarotNl3(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl3).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl3 = sym::IkResidualFuncCost3WrtPaNl3(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl3).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl3 = sym::IkResidualFuncCost3WrtPbNl3(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl3).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl3 = sym::IkResidualFuncCost3WrtPcNl3(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl3).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl3 = sym::IkResidualFuncCost3WrtPdNl3(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl3).finished();
                    }
                    break;



                case 4:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl4 = sym::IkResidualFuncCost3WrtFh1Nl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl4).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl4 = sym::IkResidualFuncCost3WrtFv1Nl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl4).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl4 = sym::IkResidualFuncCost3WrtDeltarotNl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl4).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl4 = sym::IkResidualFuncCost3WrtPaNl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl4).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl4 = sym::IkResidualFuncCost3WrtPbNl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl4).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl4 = sym::IkResidualFuncCost3WrtPcNl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl4).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl4 = sym::IkResidualFuncCost3WrtPdNl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl4).finished();
                    }
                    break;


                case 5:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl5(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl5 = sym::IkResidualFuncCost3WrtFh1Nl5(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl5).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl5 = sym::IkResidualFuncCost3WrtFv1Nl5(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl5).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl5 = sym::IkResidualFuncCost3WrtDeltarotNl5(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl5).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl5 = sym::IkResidualFuncCost3WrtPaNl5(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl5).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl5 = sym::IkResidualFuncCost3WrtPbNl5(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl5).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl5 = sym::IkResidualFuncCost3WrtPcNl5(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl5).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl5 = sym::IkResidualFuncCost3WrtPdNl5(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl5).finished();
                    }
                    break;


                case 6:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl6 = sym::IkResidualFuncCost3WrtFh1Nl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl6).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl6 = sym::IkResidualFuncCost3WrtFv1Nl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl6).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl6 = sym::IkResidualFuncCost3WrtDeltarotNl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl6).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl6 = sym::IkResidualFuncCost3WrtPaNl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl6).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl6 = sym::IkResidualFuncCost3WrtPbNl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl6).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl6 = sym::IkResidualFuncCost3WrtPcNl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl6).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl6 = sym::IkResidualFuncCost3WrtPdNl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl6).finished();
                    }
                    break;



                case 7:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl7(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl7 = sym::IkResidualFuncCost3WrtFh1Nl7(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl7).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl7 = sym::IkResidualFuncCost3WrtFv1Nl7(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl7).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl7 = sym::IkResidualFuncCost3WrtDeltarotNl7(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl7).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl7 = sym::IkResidualFuncCost3WrtPaNl7(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl7).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl7 = sym::IkResidualFuncCost3WrtPbNl7(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl7).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl7 = sym::IkResidualFuncCost3WrtPcNl7(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl7).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl7 = sym::IkResidualFuncCost3WrtPdNl7(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl7).finished();
                    }
                    break;


                case 8:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl8 = sym::IkResidualFuncCost3WrtFh1Nl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl8).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl8 = sym::IkResidualFuncCost3WrtFv1Nl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl8).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl8 = sym::IkResidualFuncCost3WrtDeltarotNl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl8).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl8 = sym::IkResidualFuncCost3WrtPaNl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl8).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl8 = sym::IkResidualFuncCost3WrtPbNl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl8).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl8 = sym::IkResidualFuncCost3WrtPcNl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl8).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl8 = sym::IkResidualFuncCost3WrtPdNl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl8).finished();
                    }
                    break;



                case 9:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl9 = sym::IkResidualFuncCost3WrtFh1Nl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl9).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl9 = sym::IkResidualFuncCost3WrtFv1Nl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl9).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl9 = sym::IkResidualFuncCost3WrtDeltarotNl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl9).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl9 = sym::IkResidualFuncCost3WrtPaNl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl9).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl9 = sym::IkResidualFuncCost3WrtPbNl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl9).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl9 = sym::IkResidualFuncCost3WrtPcNl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl9).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl9 = sym::IkResidualFuncCost3WrtPdNl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl9).finished();
                    }
                    break;



                case 10:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl10 = sym::IkResidualFuncCost3WrtFh1Nl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl10).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl10 = sym::IkResidualFuncCost3WrtFv1Nl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl10).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl10 = sym::IkResidualFuncCost3WrtDeltarotNl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl10).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl10 = sym::IkResidualFuncCost3WrtPaNl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl10).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl10 = sym::IkResidualFuncCost3WrtPbNl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl10).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl10 = sym::IkResidualFuncCost3WrtPcNl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl10).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl10 = sym::IkResidualFuncCost3WrtPdNl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl10).finished();
                    }
                    break;



                case 11:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl11 = sym::IkResidualFuncCost3WrtFh1Nl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl11).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl11 = sym::IkResidualFuncCost3WrtFv1Nl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl11).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl11 = sym::IkResidualFuncCost3WrtDeltarotNl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl11).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl11 = sym::IkResidualFuncCost3WrtPaNl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl11).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl11 = sym::IkResidualFuncCost3WrtPbNl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl11).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl11 = sym::IkResidualFuncCost3WrtPcNl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl11).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl11 = sym::IkResidualFuncCost3WrtPdNl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl11).finished();
                    }
                    break;



                case 12:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl12 = sym::IkResidualFuncCost3WrtFh1Nl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl12).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl12 = sym::IkResidualFuncCost3WrtFv1Nl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl12).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl12 = sym::IkResidualFuncCost3WrtDeltarotNl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl12).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl12 = sym::IkResidualFuncCost3WrtPaNl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl12).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl12 = sym::IkResidualFuncCost3WrtPbNl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl12).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl12 = sym::IkResidualFuncCost3WrtPcNl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl12).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl12 = sym::IkResidualFuncCost3WrtPdNl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl12).finished();
                    }
                    break;



                case 13:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl13 = sym::IkResidualFuncCost3WrtFh1Nl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl13).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl13 = sym::IkResidualFuncCost3WrtFv1Nl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl13).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl13 = sym::IkResidualFuncCost3WrtDeltarotNl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl13).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl13 = sym::IkResidualFuncCost3WrtPaNl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl13).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl13 = sym::IkResidualFuncCost3WrtPbNl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl13).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl13 = sym::IkResidualFuncCost3WrtPcNl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl13).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl13 = sym::IkResidualFuncCost3WrtPdNl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl13).finished();
                    }
                    break;



                case 14:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl14 = sym::IkResidualFuncCost3WrtFh1Nl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl14).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl14 = sym::IkResidualFuncCost3WrtFv1Nl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl14).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl14 = sym::IkResidualFuncCost3WrtDeltarotNl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl14).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl14 = sym::IkResidualFuncCost3WrtPaNl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl14).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl14 = sym::IkResidualFuncCost3WrtPbNl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl14).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl14 = sym::IkResidualFuncCost3WrtPcNl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl14).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl14 = sym::IkResidualFuncCost3WrtPdNl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl14).finished();
                    }
                    break;



                case 15:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl15 = sym::IkResidualFuncCost3WrtFh1Nl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl15).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl15 = sym::IkResidualFuncCost3WrtFv1Nl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl15).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl15 = sym::IkResidualFuncCost3WrtDeltarotNl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl15).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl15 = sym::IkResidualFuncCost3WrtPaNl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl15).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl15 = sym::IkResidualFuncCost3WrtPbNl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl15).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl15 = sym::IkResidualFuncCost3WrtPcNl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl15).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl15 = sym::IkResidualFuncCost3WrtPdNl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl15).finished();
                    }
                    break;



                case 16:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl16 = sym::IkResidualFuncCost3WrtFh1Nl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl16).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl16 = sym::IkResidualFuncCost3WrtFv1Nl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl16).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl16 = sym::IkResidualFuncCost3WrtDeltarotNl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl16).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl16 = sym::IkResidualFuncCost3WrtPaNl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl16).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl16 = sym::IkResidualFuncCost3WrtPbNl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl16).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl16 = sym::IkResidualFuncCost3WrtPcNl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl16).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl16 = sym::IkResidualFuncCost3WrtPdNl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl16).finished();
                    }
                    break;



                case 17:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl17 = sym::IkResidualFuncCost3WrtFh1Nl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl17).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl17 = sym::IkResidualFuncCost3WrtFv1Nl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl17).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl17 = sym::IkResidualFuncCost3WrtDeltarotNl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl17).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl17 = sym::IkResidualFuncCost3WrtPaNl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl17).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl17 = sym::IkResidualFuncCost3WrtPbNl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl17).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl17 = sym::IkResidualFuncCost3WrtPcNl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl17).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl17 = sym::IkResidualFuncCost3WrtPdNl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl17).finished();
                    }
                    break;


                case 18:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl18 = sym::IkResidualFuncCost3WrtFh1Nl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl18).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl18 = sym::IkResidualFuncCost3WrtFv1Nl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl18).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl18 = sym::IkResidualFuncCost3WrtDeltarotNl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl18).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl18 = sym::IkResidualFuncCost3WrtPaNl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl18).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl18 = sym::IkResidualFuncCost3WrtPbNl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl18).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl18 = sym::IkResidualFuncCost3WrtPcNl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl18).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl18 = sym::IkResidualFuncCost3WrtPdNl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl18).finished();
                    }
                    break;



                case 19:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl19 = sym::IkResidualFuncCost3WrtFh1Nl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl19).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl19 = sym::IkResidualFuncCost3WrtFv1Nl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl19).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl19 = sym::IkResidualFuncCost3WrtDeltarotNl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl19).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl19 = sym::IkResidualFuncCost3WrtPaNl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl19).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl19 = sym::IkResidualFuncCost3WrtPbNl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl19).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl19 = sym::IkResidualFuncCost3WrtPcNl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl19).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl19 = sym::IkResidualFuncCost3WrtPdNl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl19).finished();
                    }
                    break;



                case 20:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl20 = sym::IkResidualFuncCost3WrtFh1Nl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl20).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl20 = sym::IkResidualFuncCost3WrtFv1Nl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl20).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl20 = sym::IkResidualFuncCost3WrtDeltarotNl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl20).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl20 = sym::IkResidualFuncCost3WrtPaNl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl20).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl20 = sym::IkResidualFuncCost3WrtPbNl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl20).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl20 = sym::IkResidualFuncCost3WrtPcNl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl20).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl20 = sym::IkResidualFuncCost3WrtPdNl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl20).finished();
                    }
                    break;



                case 21:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl21 = sym::IkResidualFuncCost3WrtFh1Nl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl21).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl21 = sym::IkResidualFuncCost3WrtFv1Nl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl21).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl21 = sym::IkResidualFuncCost3WrtDeltarotNl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl21).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl21 = sym::IkResidualFuncCost3WrtPaNl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl21).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl21 = sym::IkResidualFuncCost3WrtPbNl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl21).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl21 = sym::IkResidualFuncCost3WrtPcNl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl21).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl21 = sym::IkResidualFuncCost3WrtPdNl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl21).finished();
                    }
                    break;



                case 22:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl22 = sym::IkResidualFuncCost3WrtFh1Nl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl22).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl22 = sym::IkResidualFuncCost3WrtFv1Nl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl22).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl22 = sym::IkResidualFuncCost3WrtDeltarotNl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl22).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl22 = sym::IkResidualFuncCost3WrtPaNl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl22).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl22 = sym::IkResidualFuncCost3WrtPbNl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl22).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl22 = sym::IkResidualFuncCost3WrtPcNl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl22).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl22 = sym::IkResidualFuncCost3WrtPdNl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl22).finished();
                    }
                    break;



                case 23:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl23 = sym::IkResidualFuncCost3WrtFh1Nl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl23).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl23 = sym::IkResidualFuncCost3WrtFv1Nl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl23).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl23 = sym::IkResidualFuncCost3WrtDeltarotNl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl23).finished();
                    }
                    if(H4)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pa_Nl23 = sym::IkResidualFuncCost3WrtPaNl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H4 = (Matrix(4, 3) << IKresidual_func_wrt_pa_Nl23).finished();
                    }
                    if(H5)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pb_Nl23 = sym::IkResidualFuncCost3WrtPbNl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H5 = (Matrix(4, 3) << IKresidual_func_wrt_pb_Nl23).finished();
                    }
                    if(H6)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pc_Nl23 = sym::IkResidualFuncCost3WrtPcNl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H6 = (Matrix(4, 3) << IKresidual_func_wrt_pc_Nl23).finished();
                    }
                    if(H7)
                    {
                        Eigen::Matrix<double, 4, 3> IKresidual_func_wrt_pd_Nl23 = sym::IkResidualFuncCost3WrtPdNl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init, p_a, p_b, p_c, p_d, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), sym::kDefaultEpsilon<double>);
                        *H7 = (Matrix(4, 3) << IKresidual_func_wrt_pd_Nl23).finished();
                    }
                    break;
            }
            // std::cout << "cost3: " << Ikresidual_func << std::endl;
            return (Vector(4) << Ikresidual_func).finished();
        }
    };
}