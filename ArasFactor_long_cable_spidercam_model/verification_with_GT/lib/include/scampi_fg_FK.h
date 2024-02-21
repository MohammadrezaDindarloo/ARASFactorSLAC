#pragma once

#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <boost/optional.hpp>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/inference/Symbol.h>
#include "scampi_function_header_include.h"

using namespace gtsam;
using namespace std;

namespace gtsam
{
    class FK_factor_graoh_cost1 : public NoiseModelFactorN<double, double, gtsam::Rot3, gtsam::Point3>
    {

    private:
        double lc0;
        double lc1; 
        double lc2; 
        double lc3;         
        double rot_init_x; 
        double rot_init_y; 
        double rot_init_z; 
        double rot_init_w; 
        double epsilon = 0.0;

    public: 
        // Constructor
        FK_factor_graoh_cost1(Key key1, Key key2, Key key3, Key key4, double lc0_, double lc1_, double lc2_, double lc3_, double rot_init_x_, double rot_init_y_, double rot_init_z_, double rot_init_w_, const SharedNoiseModel &model) 
        : NoiseModelFactorN<double, double, gtsam::Rot3, gtsam::Point3>(model, key1, key2, key3, key4), lc0(lc0_), lc1(lc1_), lc2(lc2_), lc3(lc3_), rot_init_x(rot_init_x_), rot_init_y(rot_init_y_), rot_init_z(rot_init_z_), rot_init_w(rot_init_w_) {}

        // Evaluate the error
        Vector evaluateError(const double &fh1, const double &fv1, const gtsam::Rot3 &DeltaRot, const gtsam::Point3 &position_vector,
                             OptionalMatrixType H1,
                             OptionalMatrixType H2,
                             OptionalMatrixType H3,
                             OptionalMatrixType H4) const override
        {   
            Eigen::Matrix<double, 4, 1> Fkresidual_func = sym::FkResidualFuncCost1(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, lc0, lc1, lc2, lc3, epsilon);
            if(H1)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fh1 = sym::FkResidualFuncCost1WrtFh1(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, lc0, lc1, lc2, lc3, epsilon);
                *H1 = (Matrix(4, 1) << Fkresidual_func_wrt_fh1).finished();
            }
            if(H2)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fv1 = sym::FkResidualFuncCost1WrtFv1(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, lc0, lc1, lc2, lc3, epsilon);
                *H2 = (Matrix(4, 1) << Fkresidual_func_wrt_fv1).finished();
            }
            if(H3) 
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_Deltarot = sym::FkResidualFuncCost1WrtDeltarot(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, lc0, lc1, lc2, lc3, epsilon);
                *H3 = (Matrix(4, 3) << Fkresidual_func_wrt_Deltarot).finished();
            }
            if(H4)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_PositionVector = sym::FkResidualFuncCost1WrtPositionVector(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, lc0, lc1, lc2, lc3, epsilon);
                *H4 = (Matrix(4, 3) << Fkresidual_func_wrt_PositionVector).finished();
            }

            return (Vector(4) << Fkresidual_func).finished();
        }
    };


    class FK_factor_graoh_cost2 : public NoiseModelFactorN<double, double, gtsam::Rot3, gtsam::Point3>
    {

    private:
        double lc0;
        double lc1; 
        double lc2; 
        double lc3;         
        double rot_init_x; 
        double rot_init_y; 
        double rot_init_z; 
        double rot_init_w; 
        double epsilon = 0.0;

    public:
        // Constructor
        FK_factor_graoh_cost2(Key key1, Key key2, Key key3, Key key4, double lc0_, double lc1_, double lc2_, double lc3_, double rot_init_x_, double rot_init_y_, double rot_init_z_, double rot_init_w_, const SharedNoiseModel &model) 
        : NoiseModelFactorN<double, double, gtsam::Rot3, gtsam::Point3>(model, key1, key2, key3, key4), lc0(lc0_), lc1(lc1_), lc2(lc2_), lc3(lc3_), rot_init_x(rot_init_x_), rot_init_y(rot_init_y_), rot_init_z(rot_init_z_), rot_init_w(rot_init_w_) {}

        // Evaluate the error
        Vector evaluateError(const double &fh1, const double &fv1, const gtsam::Rot3 &DeltaRot, const gtsam::Point3 &position_vector,
                             OptionalMatrixType H1,
                             OptionalMatrixType H2,
                             OptionalMatrixType H3,
                             OptionalMatrixType H4) const override
        {   
            Eigen::Matrix<double, 4, 1> Fkresidual_func = sym::FkResidualFuncCost2(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, lc0, lc1, lc2, lc3, epsilon);
            if(H1)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fh1 = sym::FkResidualFuncCost2WrtFh1(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, lc0, lc1, lc2, lc3, epsilon);
                *H1 = (Matrix(4, 1) << Fkresidual_func_wrt_fh1).finished();
            }
            if(H2)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fv1 = sym::FkResidualFuncCost2WrtFv1(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, lc0, lc1, lc2, lc3, epsilon);
                *H2 = (Matrix(4, 1) << Fkresidual_func_wrt_fv1).finished();
            }
            if(H3)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_Deltarot = sym::FkResidualFuncCost2WrtDeltarot(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, lc0, lc1, lc2, lc3, epsilon);
                *H3 = (Matrix(4, 3) << Fkresidual_func_wrt_Deltarot).finished();
            }
            if(H4)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_PositionVector = sym::FkResidualFuncCost2WrtPositionVector(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, lc0, lc1, lc2, lc3, epsilon);
                *H4 = (Matrix(4, 3) << Fkresidual_func_wrt_PositionVector).finished();
            }

            return (Vector(4) << Fkresidual_func).finished();
        }
    };


    class FK_factor_graoh_cost3 : public NoiseModelFactorN<double, double, gtsam::Rot3, gtsam::Point3>
    {

    private:
        double lc0;
        double lc1; 
        double lc2; 
        double lc3;         
        double rot_init_x; 
        double rot_init_y; 
        double rot_init_z; 
        double rot_init_w; 
        double epsilon = 0.0;

    public:
        // Constructor
        FK_factor_graoh_cost3(Key key1, Key key2, Key key3, Key key4, double lc0_, double lc1_, double lc2_, double lc3_, double rot_init_x_, double rot_init_y_, double rot_init_z_, double rot_init_w_, const SharedNoiseModel &model) 
        : NoiseModelFactorN<double, double, gtsam::Rot3, gtsam::Point3>(model, key1, key2, key3, key4), lc0(lc0_), lc1(lc1_), lc2(lc2_), lc3(lc3_), rot_init_x(rot_init_x_), rot_init_y(rot_init_y_), rot_init_z(rot_init_z_), rot_init_w(rot_init_w_) {}

        // Evaluate the error
        Vector evaluateError(const double &fh1, const double &fv1, const gtsam::Rot3 &DeltaRot, const gtsam::Point3 &position_vector,
                             OptionalMatrixType H1,
                             OptionalMatrixType H2,
                             OptionalMatrixType H3,
                             OptionalMatrixType H4) const override
        {   
            Eigen::Matrix<double, 4, 1> Fkresidual_func = sym::FkResidualFuncCost3(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, lc0, lc1, lc2, lc3, epsilon);
            if(H1)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fh1 = sym::FkResidualFuncCost3WrtFh1(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, lc0, lc1, lc2, lc3, epsilon);
                *H1 = (Matrix(4, 1) << Fkresidual_func_wrt_fh1).finished();
            }
            if(H2)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fv1 = sym::FkResidualFuncCost3WrtFv1(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, lc0, lc1, lc2, lc3, epsilon);
                *H2 = (Matrix(4, 1) << Fkresidual_func_wrt_fv1).finished();
            }
            if(H3)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_Deltarot = sym::FkResidualFuncCost3WrtDeltarot(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, lc0, lc1, lc2, lc3, epsilon);
                *H3 = (Matrix(4, 3) << Fkresidual_func_wrt_Deltarot).finished();
            }
            if(H4)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_PositionVector = sym::FkResidualFuncCost3WrtPositionVector(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, lc0, lc1, lc2, lc3, epsilon);
                *H4 = (Matrix(4, 3) << Fkresidual_func_wrt_PositionVector).finished();
            }

            return (Vector(4) << Fkresidual_func).finished();
        }
    };

}
