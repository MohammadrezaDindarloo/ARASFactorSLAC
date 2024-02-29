#pragma once

#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <boost/optional.hpp>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/inference/Symbol.h>
#include "scampi_function_header_include.h"
#include <gtsam/base/Matrix.h>
#include <gtsam/geometry/Pose3.h>
#include <gtsam/geometry/Rot3.h>
#include <sym/rot3.h>


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

namespace gtsam
{
    class IK_factor_graoh_cost1 : public NoiseModelFactorN<double, double, gtsam::Rot3>
    {

    private:
        double p_init0;
        double p_init1; 
        double p_init2; 
        double rot_init_x; 
        double rot_init_y; 
        double rot_init_z; 
        double rot_init_w; 
        double epsilon = 0.0;
        int largest_cable = 0;

    public:
        // Constructor
        IK_factor_graoh_cost1(Key key1, Key key2, Key key3, double p_init0_, double p_init1_, double p_init2_, double rot_init_x_, double rot_init_y_, double rot_init_z_, double rot_init_w_, const int largest_cable_, const SharedNoiseModel &model) 
        : NoiseModelFactorN<double, double, gtsam::Rot3>(model, key1, key2, key3), p_init0(p_init0_), p_init1(p_init1_), p_init2(p_init2_), rot_init_x(rot_init_x_), rot_init_y(rot_init_y_), rot_init_z(rot_init_z_), rot_init_w(rot_init_w_), largest_cable(largest_cable_) {}
        
        // Evaluate the error
        Vector evaluateError(const double &fh1, const double &fv1, const gtsam::Rot3 &DeltaRot,
                             OptionalMatrixType H1,
                             OptionalMatrixType H2,
                             OptionalMatrixType H3) const override
        {   
            Eigen::Matrix<double, 4, 1> Ikresidual_func;

            switch (largest_cable)

            {   

                case 0:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl0 = sym::IkResidualFuncCost1WrtFh1Nl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl0).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl0 = sym::IkResidualFuncCost1WrtFv1Nl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl0).finished();
                    }
                    if (H3) 
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl0 = sym::IkResidualFuncCost1WrtDeltarotNl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl0).finished();
                    }
                    break;


                case 1:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl1 = sym::IkResidualFuncCost1WrtFh1Nl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl1).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl1 = sym::IkResidualFuncCost1WrtFv1Nl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl1).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl1 = sym::IkResidualFuncCost1WrtDeltarotNl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl1).finished();
                    }
                    break;


                case 2:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl2 = sym::IkResidualFuncCost1WrtFh1Nl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl2).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl2 = sym::IkResidualFuncCost1WrtFv1Nl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl2).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl2 = sym::IkResidualFuncCost1WrtDeltarotNl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl2).finished();
                    }
                    break;



                case 3:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl1 = sym::IkResidualFuncCost1WrtFh1Nl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl1).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl1 = sym::IkResidualFuncCost1WrtFv1Nl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl1).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl1 = sym::IkResidualFuncCost1WrtDeltarotNl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl1).finished();
                    }
                    break;



                case 4:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl4 = sym::IkResidualFuncCost1WrtFh1Nl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl4).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl4 = sym::IkResidualFuncCost1WrtFv1Nl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl4).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl4 = sym::IkResidualFuncCost1WrtDeltarotNl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl4).finished();
                    }
                    break;



                case 5:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl6 = sym::IkResidualFuncCost1WrtFh1Nl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl6).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl6 = sym::IkResidualFuncCost1WrtFv1Nl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl6).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl6 = sym::IkResidualFuncCost1WrtDeltarotNl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl6).finished();
                    }
                    break;



                case 7:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl7 = sym::IkResidualFuncCost1WrtFh1Nl7(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl7).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl7 = sym::IkResidualFuncCost1WrtFv1Nl7(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl7).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl7 = sym::IkResidualFuncCost1WrtDeltarotNl7(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl7).finished();
                    }
                    break;


                case 8:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl8 = sym::IkResidualFuncCost1WrtFh1Nl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl8).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl8 = sym::IkResidualFuncCost1WrtFv1Nl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl8).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl8 = sym::IkResidualFuncCost1WrtDeltarotNl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl8).finished();
                    }
                    break;



                case 9:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl9 = sym::IkResidualFuncCost1WrtFh1Nl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl9).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl9 = sym::IkResidualFuncCost1WrtFv1Nl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl9).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl9 = sym::IkResidualFuncCost1WrtDeltarotNl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl9).finished();
                    }
                    break;



                case 10:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl10 = sym::IkResidualFuncCost1WrtFh1Nl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl10).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl10 = sym::IkResidualFuncCost1WrtFv1Nl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl10).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl10 = sym::IkResidualFuncCost1WrtDeltarotNl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl10).finished();
                    }
                    break;



                case 11:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl11 = sym::IkResidualFuncCost1WrtFh1Nl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl11).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl11 = sym::IkResidualFuncCost1WrtFv1Nl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl11).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl11 = sym::IkResidualFuncCost1WrtDeltarotNl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl11).finished();
                    }
                    break;



                case 12:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl12 = sym::IkResidualFuncCost1WrtFh1Nl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl12).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl12 = sym::IkResidualFuncCost1WrtFv1Nl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl12).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl12 = sym::IkResidualFuncCost1WrtDeltarotNl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl12).finished();
                    }
                    break;



                case 13:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl13 = sym::IkResidualFuncCost1WrtFh1Nl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl13).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl13 = sym::IkResidualFuncCost1WrtFv1Nl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl13).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl13 = sym::IkResidualFuncCost1WrtDeltarotNl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl13).finished();
                    }
                    break;



                case 14:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl14 = sym::IkResidualFuncCost1WrtFh1Nl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl14).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl14 = sym::IkResidualFuncCost1WrtFv1Nl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl14).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl14 = sym::IkResidualFuncCost1WrtDeltarotNl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl14).finished();
                    }
                    break;



                case 15:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl15 = sym::IkResidualFuncCost1WrtFh1Nl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl15).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl15 = sym::IkResidualFuncCost1WrtFv1Nl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl15).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl15 = sym::IkResidualFuncCost1WrtDeltarotNl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl15).finished();
                    }
                    break;



                case 16:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl16 = sym::IkResidualFuncCost1WrtFh1Nl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl16).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl16 = sym::IkResidualFuncCost1WrtFv1Nl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl16).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl16 = sym::IkResidualFuncCost1WrtDeltarotNl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl16).finished();
                    }
                    break;



                case 17:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl17 = sym::IkResidualFuncCost1WrtFh1Nl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl17).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl17 = sym::IkResidualFuncCost1WrtFv1Nl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl17).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl17 = sym::IkResidualFuncCost1WrtDeltarotNl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl17).finished();
                    }
                    break;


                case 18:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl18 = sym::IkResidualFuncCost1WrtFh1Nl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl18).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl18 = sym::IkResidualFuncCost1WrtFv1Nl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl18).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl18 = sym::IkResidualFuncCost1WrtDeltarotNl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl18).finished();
                    }
                    break;



                case 19:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl19 = sym::IkResidualFuncCost1WrtFh1Nl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl19).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl19 = sym::IkResidualFuncCost1WrtFv1Nl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl19).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl19 = sym::IkResidualFuncCost1WrtDeltarotNl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl19).finished();
                    }
                    break;



                case 20:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl20 = sym::IkResidualFuncCost1WrtFh1Nl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl20).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl20 = sym::IkResidualFuncCost1WrtFv1Nl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl20).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl20 = sym::IkResidualFuncCost1WrtDeltarotNl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl20).finished();
                    }
                    break;



                case 21:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl21 = sym::IkResidualFuncCost1WrtFh1Nl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl21).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl21 = sym::IkResidualFuncCost1WrtFv1Nl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl21).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl21 = sym::IkResidualFuncCost1WrtDeltarotNl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl21).finished();
                    }
                    break;



                case 22:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl22 = sym::IkResidualFuncCost1WrtFh1Nl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl22).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl22 = sym::IkResidualFuncCost1WrtFv1Nl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl22).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl22 = sym::IkResidualFuncCost1WrtDeltarotNl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl22).finished();
                    }
                    break;



                case 23:
                    Ikresidual_func = sym::IkResidualFuncCost1Nl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl23 = sym::IkResidualFuncCost1WrtFh1Nl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl23).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl23 = sym::IkResidualFuncCost1WrtFv1Nl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl23).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl23 = sym::IkResidualFuncCost1WrtDeltarotNl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl23).finished();
                    }
                    break;
            }
            return (Vector(4) << Ikresidual_func).finished();
        }
    };

    class IK_factor_graoh_cost2 : public NoiseModelFactorN<double, double, gtsam::Rot3>
    {

    private:
        double p_init0;
        double p_init1; 
        double p_init2; 
        double rot_init_x; 
        double rot_init_y; 
        double rot_init_z; 
        double rot_init_w; 
        double epsilon = 0.0;
        int largest_cable = 0;

    public:
        // Constructor
        IK_factor_graoh_cost2(Key key1, Key key2, Key key3, double p_init0_, double p_init1_, double p_init2_, double rot_init_x_, double rot_init_y_, double rot_init_z_, double rot_init_w_, const int largest_cable_, const SharedNoiseModel &model) 
        : NoiseModelFactorN<double, double, gtsam::Rot3>(model, key1, key2, key3), p_init0(p_init0_), p_init1(p_init1_), p_init2(p_init2_), rot_init_x(rot_init_x_), rot_init_y(rot_init_y_), rot_init_z(rot_init_z_), rot_init_w(rot_init_w_), largest_cable(largest_cable_) {}

        // Evaluate the error
        Vector evaluateError(const double &fh1, const double &fv1, const gtsam::Rot3 &DeltaRot,
                             OptionalMatrixType H1,
                             OptionalMatrixType H2,
                             OptionalMatrixType H3) const override
        {   
            Eigen::Matrix<double, 4, 1> Ikresidual_func;

            switch (largest_cable)
            {   
                case 0:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl0 = sym::IkResidualFuncCost2WrtFh1Nl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl0).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl0 = sym::IkResidualFuncCost2WrtFv1Nl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl0).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl0 = sym::IkResidualFuncCost2WrtDeltarotNl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl0).finished();
                    }
                    break;


                case 1:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl1 = sym::IkResidualFuncCost2WrtFh1Nl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl1).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl1 = sym::IkResidualFuncCost2WrtFv1Nl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl1).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl1 = sym::IkResidualFuncCost2WrtDeltarotNl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl1).finished();
                    }
                    break;


                case 2:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl2 = sym::IkResidualFuncCost2WrtFh1Nl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl2).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl2 = sym::IkResidualFuncCost2WrtFv1Nl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl2).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl2 = sym::IkResidualFuncCost2WrtDeltarotNl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl2).finished();
                    }
                    break;



                case 3:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl1 = sym::IkResidualFuncCost2WrtFh1Nl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl1).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl1 = sym::IkResidualFuncCost2WrtFv1Nl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl1).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl1 = sym::IkResidualFuncCost2WrtDeltarotNl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl1).finished();
                    }
                    break;



                case 4:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl4 = sym::IkResidualFuncCost2WrtFh1Nl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl4).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl4 = sym::IkResidualFuncCost2WrtFv1Nl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl4).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl4 = sym::IkResidualFuncCost2WrtDeltarotNl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl4).finished();
                    }
                    break;



                case 5:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl6 = sym::IkResidualFuncCost2WrtFh1Nl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl6).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl6 = sym::IkResidualFuncCost2WrtFv1Nl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl6).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl6 = sym::IkResidualFuncCost2WrtDeltarotNl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl6).finished();
                    }
                    break;



                case 7:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl7 = sym::IkResidualFuncCost2WrtFh1Nl7(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl7).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl7 = sym::IkResidualFuncCost2WrtFv1Nl7(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl7).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl7 = sym::IkResidualFuncCost2WrtDeltarotNl7(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl7).finished();
                    }
                    break;


                case 8:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl8 = sym::IkResidualFuncCost2WrtFh1Nl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl8).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl8 = sym::IkResidualFuncCost2WrtFv1Nl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl8).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl8 = sym::IkResidualFuncCost2WrtDeltarotNl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl8).finished();
                    }
                    break;



                case 9:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl9 = sym::IkResidualFuncCost2WrtFh1Nl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl9).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl9 = sym::IkResidualFuncCost2WrtFv1Nl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl9).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl9 = sym::IkResidualFuncCost2WrtDeltarotNl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl9).finished();
                    }
                    break;



                case 10:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl10 = sym::IkResidualFuncCost2WrtFh1Nl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl10).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl10 = sym::IkResidualFuncCost2WrtFv1Nl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl10).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl10 = sym::IkResidualFuncCost2WrtDeltarotNl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl10).finished();
                    }
                    break;



                case 11:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl11 = sym::IkResidualFuncCost2WrtFh1Nl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl11).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl11 = sym::IkResidualFuncCost2WrtFv1Nl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl11).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl11 = sym::IkResidualFuncCost2WrtDeltarotNl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl11).finished();
                    }
                    break;



                case 12:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl12 = sym::IkResidualFuncCost2WrtFh1Nl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl12).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl12 = sym::IkResidualFuncCost2WrtFv1Nl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl12).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl12 = sym::IkResidualFuncCost2WrtDeltarotNl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl12).finished();
                    }
                    break;



                case 13:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl13 = sym::IkResidualFuncCost2WrtFh1Nl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl13).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl13 = sym::IkResidualFuncCost2WrtFv1Nl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl13).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl13 = sym::IkResidualFuncCost2WrtDeltarotNl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl13).finished();
                    }
                    break;



                case 14:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl14 = sym::IkResidualFuncCost2WrtFh1Nl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl14).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl14 = sym::IkResidualFuncCost2WrtFv1Nl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl14).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl14 = sym::IkResidualFuncCost2WrtDeltarotNl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl14).finished();
                    }
                    break;



                case 15:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl15 = sym::IkResidualFuncCost2WrtFh1Nl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl15).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl15 = sym::IkResidualFuncCost2WrtFv1Nl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl15).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl15 = sym::IkResidualFuncCost2WrtDeltarotNl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl15).finished();
                    }
                    break;



                case 16:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl16 = sym::IkResidualFuncCost2WrtFh1Nl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl16).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl16 = sym::IkResidualFuncCost2WrtFv1Nl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl16).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl16 = sym::IkResidualFuncCost2WrtDeltarotNl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl16).finished();
                    }
                    break;



                case 17:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl17 = sym::IkResidualFuncCost2WrtFh1Nl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl17).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl17 = sym::IkResidualFuncCost2WrtFv1Nl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl17).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl17 = sym::IkResidualFuncCost2WrtDeltarotNl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl17).finished();
                    }
                    break;


                case 18:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl18 = sym::IkResidualFuncCost2WrtFh1Nl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl18).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl18 = sym::IkResidualFuncCost2WrtFv1Nl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl18).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl18 = sym::IkResidualFuncCost2WrtDeltarotNl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl18).finished();
                    }
                    break;



                case 19:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl19 = sym::IkResidualFuncCost2WrtFh1Nl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl19).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl19 = sym::IkResidualFuncCost2WrtFv1Nl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl19).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl19 = sym::IkResidualFuncCost2WrtDeltarotNl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl19).finished();
                    }
                    break;



                case 20:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl20 = sym::IkResidualFuncCost2WrtFh1Nl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl20).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl20 = sym::IkResidualFuncCost2WrtFv1Nl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl20).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl20 = sym::IkResidualFuncCost2WrtDeltarotNl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl20).finished();
                    }
                    break;



                case 21:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl21 = sym::IkResidualFuncCost2WrtFh1Nl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl21).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl21 = sym::IkResidualFuncCost2WrtFv1Nl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl21).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl21 = sym::IkResidualFuncCost2WrtDeltarotNl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl21).finished();
                    }
                    break;



                case 22:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl22 = sym::IkResidualFuncCost2WrtFh1Nl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl22).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl22 = sym::IkResidualFuncCost2WrtFv1Nl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl22).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl22 = sym::IkResidualFuncCost2WrtDeltarotNl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl22).finished();
                    }
                    break;



                case 23:
                    Ikresidual_func = sym::IkResidualFuncCost2Nl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl23 = sym::IkResidualFuncCost2WrtFh1Nl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl23).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl23 = sym::IkResidualFuncCost2WrtFv1Nl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl23).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl23 = sym::IkResidualFuncCost2WrtDeltarotNl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl23).finished();
                    }
                    break;
            }
            return (Vector(4) << Ikresidual_func).finished(); 
        }
    };


    class IK_factor_graoh_cost3 : public NoiseModelFactorN<double, double, gtsam::Rot3>
    {

    private:
        double p_init0;
        double p_init1; 
        double p_init2; 
        double rot_init_x; 
        double rot_init_y; 
        double rot_init_z; 
        double rot_init_w; 
        double epsilon = 0.0;
        int largest_cable;

    public:
        // Constructor
        IK_factor_graoh_cost3(Key key1, Key key2, Key key3, double p_init0_, double p_init1_, double p_init2_, double rot_init_x_, double rot_init_y_, double rot_init_z_, double rot_init_w_, const int largest_cable_, const SharedNoiseModel &model) 
        : NoiseModelFactorN<double, double, gtsam::Rot3>(model, key1, key2, key3), p_init0(p_init0_), p_init1(p_init1_), p_init2(p_init2_), rot_init_x(rot_init_x_), rot_init_y(rot_init_y_), rot_init_z(rot_init_z_), rot_init_w(rot_init_w_), largest_cable(largest_cable_) {}

        // Evaluate the error
        Vector evaluateError(const double &fh1, const double &fv1, const gtsam::Rot3 &DeltaRot,
                             OptionalMatrixType H1,
                             OptionalMatrixType H2,
                             OptionalMatrixType H3) const override
        {   
            Eigen::Matrix<double, 4, 1> Ikresidual_func; 

            switch (largest_cable)
            {   
                case 0:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl0 = sym::IkResidualFuncCost3WrtFh1Nl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl0).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl0 = sym::IkResidualFuncCost3WrtFv1Nl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl0).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl0 = sym::IkResidualFuncCost3WrtDeltarotNl0(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl0).finished();
                    }
                    break;


                case 1:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl1 = sym::IkResidualFuncCost3WrtFh1Nl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl1).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl1 = sym::IkResidualFuncCost3WrtFv1Nl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl1).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl1 = sym::IkResidualFuncCost3WrtDeltarotNl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl1).finished();
                    }
                    break;


                case 2:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl2 = sym::IkResidualFuncCost3WrtFh1Nl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl2).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl2 = sym::IkResidualFuncCost3WrtFv1Nl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl2).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl2 = sym::IkResidualFuncCost3WrtDeltarotNl2(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl2).finished();
                    }
                    break;



                case 3:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl1 = sym::IkResidualFuncCost3WrtFh1Nl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl1).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl1 = sym::IkResidualFuncCost3WrtFv1Nl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl1).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl1 = sym::IkResidualFuncCost3WrtDeltarotNl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl1).finished();
                    }
                    break;



                case 4:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl4 = sym::IkResidualFuncCost3WrtFh1Nl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl4).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl4 = sym::IkResidualFuncCost3WrtFv1Nl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl4).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl4 = sym::IkResidualFuncCost3WrtDeltarotNl4(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl4).finished();
                    }
                    break;



                case 5:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl6 = sym::IkResidualFuncCost3WrtFh1Nl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl6).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl6 = sym::IkResidualFuncCost3WrtFv1Nl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl6).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl6 = sym::IkResidualFuncCost3WrtDeltarotNl6(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl6).finished();
                    }
                    break;



                case 7:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl1(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl7 = sym::IkResidualFuncCost3WrtFh1Nl7(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl7).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl7 = sym::IkResidualFuncCost3WrtFv1Nl7(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl7).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl7 = sym::IkResidualFuncCost3WrtDeltarotNl7(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl7).finished();
                    }
                    break;


                case 8:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl8 = sym::IkResidualFuncCost3WrtFh1Nl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl8).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl8 = sym::IkResidualFuncCost3WrtFv1Nl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl8).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl8 = sym::IkResidualFuncCost3WrtDeltarotNl8(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl8).finished();
                    }
                    break;



                case 9:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl9 = sym::IkResidualFuncCost3WrtFh1Nl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl9).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl9 = sym::IkResidualFuncCost3WrtFv1Nl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl9).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl9 = sym::IkResidualFuncCost3WrtDeltarotNl9(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl9).finished();
                    }
                    break;



                case 10:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl10 = sym::IkResidualFuncCost3WrtFh1Nl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl10).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl10 = sym::IkResidualFuncCost3WrtFv1Nl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl10).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl10 = sym::IkResidualFuncCost3WrtDeltarotNl10(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl10).finished();
                    }
                    break;



                case 11:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl11 = sym::IkResidualFuncCost3WrtFh1Nl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl11).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl11 = sym::IkResidualFuncCost3WrtFv1Nl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl11).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl11 = sym::IkResidualFuncCost3WrtDeltarotNl11(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl11).finished();
                    }
                    break;



                case 12:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl12 = sym::IkResidualFuncCost3WrtFh1Nl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl12).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl12 = sym::IkResidualFuncCost3WrtFv1Nl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl12).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl12 = sym::IkResidualFuncCost3WrtDeltarotNl12(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl12).finished();
                    }
                    break;



                case 13:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl13 = sym::IkResidualFuncCost3WrtFh1Nl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl13).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl13 = sym::IkResidualFuncCost3WrtFv1Nl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl13).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl13 = sym::IkResidualFuncCost3WrtDeltarotNl13(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl13).finished();
                    }
                    break;



                case 14:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl14 = sym::IkResidualFuncCost3WrtFh1Nl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl14).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl14 = sym::IkResidualFuncCost3WrtFv1Nl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl14).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl14 = sym::IkResidualFuncCost3WrtDeltarotNl14(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl14).finished();
                    }
                    break;



                case 15:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl15 = sym::IkResidualFuncCost3WrtFh1Nl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl15).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl15 = sym::IkResidualFuncCost3WrtFv1Nl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl15).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl15 = sym::IkResidualFuncCost3WrtDeltarotNl15(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl15).finished();
                    }
                    break;



                case 16:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl16 = sym::IkResidualFuncCost3WrtFh1Nl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl16).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl16 = sym::IkResidualFuncCost3WrtFv1Nl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl16).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl16 = sym::IkResidualFuncCost3WrtDeltarotNl16(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl16).finished();
                    }
                    break;



                case 17:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl17 = sym::IkResidualFuncCost3WrtFh1Nl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl17).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl17 = sym::IkResidualFuncCost3WrtFv1Nl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl17).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl17 = sym::IkResidualFuncCost3WrtDeltarotNl17(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl17).finished();
                    }
                    break;


                case 18:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl18 = sym::IkResidualFuncCost3WrtFh1Nl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl18).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl18 = sym::IkResidualFuncCost3WrtFv1Nl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl18).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl18 = sym::IkResidualFuncCost3WrtDeltarotNl18(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl18).finished();
                    }
                    break;



                case 19:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl19 = sym::IkResidualFuncCost3WrtFh1Nl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl19).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl19 = sym::IkResidualFuncCost3WrtFv1Nl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl19).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl19 = sym::IkResidualFuncCost3WrtDeltarotNl19(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl19).finished();
                    }
                    break;



                case 20:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl20 = sym::IkResidualFuncCost3WrtFh1Nl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl20).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl20 = sym::IkResidualFuncCost3WrtFv1Nl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl20).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl20 = sym::IkResidualFuncCost3WrtDeltarotNl20(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl20).finished();
                    }
                    break;



                case 21:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl21 = sym::IkResidualFuncCost3WrtFh1Nl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl21).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl21 = sym::IkResidualFuncCost3WrtFv1Nl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl21).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl21 = sym::IkResidualFuncCost3WrtDeltarotNl21(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl21).finished();
                    }
                    break;



                case 22:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl22 = sym::IkResidualFuncCost3WrtFh1Nl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl22).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl22 = sym::IkResidualFuncCost3WrtFv1Nl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl22).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl22 = sym::IkResidualFuncCost3WrtDeltarotNl22(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl22).finished();
                    }
                    break;



                case 23:
                    Ikresidual_func = sym::IkResidualFuncCost3Nl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                    if (H1)
                    {
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fh1_Nl23 = sym::IkResidualFuncCost3WrtFh1Nl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H1 = (Matrix(4, 1) << Ikresidual_func_wrt_fh1_Nl23).finished();
                    }
                    if (H2)
                    {   
                        Eigen::Matrix<double, 4, 1> Ikresidual_func_wrt_fv1_Nl23 = sym::IkResidualFuncCost3WrtFv1Nl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H2 = (Matrix(4, 1) << Ikresidual_func_wrt_fv1_Nl23).finished();
                    }
                    if (H3)
                    {   
                        Eigen::Matrix<double, 4, 3> Ikresidual_func_wrt_Deltarot_Nl23 = sym::IkResidualFuncCost3WrtDeltarotNl23(fh1, fv1, SymforceFromGtsam(DeltaRot), p_init0, p_init1, p_init2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, epsilon);
                        *H3 = (Matrix(4, 3) << Ikresidual_func_wrt_Deltarot_Nl23).finished();
                    }
                    break;
            }
            return (Vector(4) << Ikresidual_func).finished();
        }
    };
}