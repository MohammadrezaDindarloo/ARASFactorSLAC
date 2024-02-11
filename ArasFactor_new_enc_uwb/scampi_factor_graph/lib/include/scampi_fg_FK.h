#pragma once

#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <boost/optional.hpp>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/inference/Symbol.h>
#include "scampi_function_header_include.h"
#include "scampi_fg_IK.h"
#include <sym/rot3.h>
#include <sym/pose3.h>

using namespace gtsam;
using namespace std;

namespace gtsam
{
    class FK_factor_graoh_cost1 : public NoiseModelFactorN<double, double, gtsam::Rot3, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Point3>
    {

    private:
        Eigen::Vector4d encoder;    
        Eigen::Matrix3d rot_init;

    public: 
        // Constructor
        FK_factor_graoh_cost1(Key key1, Key key2, Key key3, Key key4, Key key5, Key key6, Key key7, Key key8, Eigen::Vector4d encoder_, Eigen::Matrix3d rot_init_, const SharedNoiseModel &model) 
        : NoiseModelFactorN<double, double, gtsam::Rot3, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Point3>(model, key1, key2, key3, key4, key5, key6, key7, key8), encoder(encoder_), rot_init(rot_init_) {}

        // Evaluate the error
        Vector evaluateError(const double &fh1, const double &fv1, const gtsam::Rot3 &DeltaRot, const gtsam::Point3 &position_vector, const gtsam::Point3 &p_a, const gtsam::Point3 &p_b, const gtsam::Point3 &p_c, const gtsam::Point3 &p_d,
                             OptionalMatrixType H1,
                             OptionalMatrixType H2,
                             OptionalMatrixType H3,
                             OptionalMatrixType H4,
                             OptionalMatrixType H5,
                             OptionalMatrixType H6,
                             OptionalMatrixType H7,
                             OptionalMatrixType H8) const override
        {   
            Eigen::Matrix<double, 4, 1> Fkresidual_func = sym::FkResidualFuncCost1(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), encoder, p_a, p_b, p_c, p_d, sym::kDefaultEpsilon<double>);
            if(H1)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fh1 = sym::FkResidualFuncCost1WrtFh1(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), encoder, p_a, p_b, p_c, p_d, sym::kDefaultEpsilon<double>);
                *H1 = (Matrix(4, 1) << Fkresidual_func_wrt_fh1).finished();
            }
            if(H2)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fv1 = sym::FkResidualFuncCost1WrtFv1(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), encoder, p_a, p_b, p_c, p_d, sym::kDefaultEpsilon<double>);
                *H2 = (Matrix(4, 1) << Fkresidual_func_wrt_fv1).finished();
            }
            if(H3) 
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_Deltarot = sym::FkResidualFuncCost1WrtDeltarot(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), encoder, p_a, p_b, p_c, p_d, sym::kDefaultEpsilon<double>);
                *H3 = (Matrix(4, 3) << Fkresidual_func_wrt_Deltarot).finished();
            }
            if(H4)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_PositionVector = sym::FkResidualFuncCost1WrtPositionVector(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), encoder, p_a, p_b, p_c, p_d, sym::kDefaultEpsilon<double>);
                *H4 = (Matrix(4, 3) << Fkresidual_func_wrt_PositionVector).finished();
            }
            if(H5)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pa = sym::FkResidualFuncCost1WrtPa(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), encoder, p_a, p_b, p_c, p_d, sym::kDefaultEpsilon<double>);
                *H5 = (Matrix(4, 3) << Fkresidual_func_wrt_pa).finished();
            }
            if(H6)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pb = sym::FkResidualFuncCost1WrtPb(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), encoder, p_a, p_b, p_c, p_d, sym::kDefaultEpsilon<double>);
                *H6 = (Matrix(4, 3) << Fkresidual_func_wrt_pb).finished();
            }
            if(H7)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pc = sym::FkResidualFuncCost1WrtPc(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), encoder, p_a, p_b, p_c, p_d, sym::kDefaultEpsilon<double>);
                *H7 = (Matrix(4, 3) << Fkresidual_func_wrt_pc).finished();
            }
            if(H8)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pd = sym::FkResidualFuncCost1WrtPd(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), encoder, p_a, p_b, p_c, p_d, sym::kDefaultEpsilon<double>);
                *H8 = (Matrix(4, 3) << Fkresidual_func_wrt_pd).finished();
            }
            // std::cout << "Fkresidual_func" << Fkresidual_func << std::endl;
            return (Vector(4) << Fkresidual_func).finished();
        }
    };


    class FK_factor_graoh_cost2 : public NoiseModelFactorN<double, double, gtsam::Rot3, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Point3>
    {

    private:
        Eigen::Vector4d encoder;    
        Eigen::Matrix3d rot_init;

    public:
        // Constructor
        FK_factor_graoh_cost2(Key key1, Key key2, Key key3, Key key4, Key key5, Key key6, Key key7, Key key8, Eigen::Vector4d encoder_, Eigen::Matrix3d rot_init_, const SharedNoiseModel &model) 
        : NoiseModelFactorN<double, double, gtsam::Rot3, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Point3>(model, key1, key2, key3, key4, key5, key6, key7, key8), encoder(encoder_), rot_init(rot_init_) {}

        // Evaluate the error
        Vector evaluateError(const double &fh1, const double &fv1, const gtsam::Rot3 &DeltaRot, const gtsam::Point3 &position_vector, const gtsam::Point3 &p_a, const gtsam::Point3 &p_b, const gtsam::Point3 &p_c, const gtsam::Point3 &p_d,
                             OptionalMatrixType H1,
                             OptionalMatrixType H2,
                             OptionalMatrixType H3,
                             OptionalMatrixType H4,
                             OptionalMatrixType H5,
                             OptionalMatrixType H6,
                             OptionalMatrixType H7,
                             OptionalMatrixType H8) const override
        {   
            Eigen::Matrix<double, 4, 1> Fkresidual_func = sym::FkResidualFuncCost2(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), encoder, p_a, p_b, p_c, p_d, sym::kDefaultEpsilon<double>);
            if(H1)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fh1 = sym::FkResidualFuncCost2WrtFh1(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), encoder, p_a, p_b, p_c, p_d, sym::kDefaultEpsilon<double>);
                *H1 = (Matrix(4, 1) << Fkresidual_func_wrt_fh1).finished();
            }
            if(H2)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fv1 = sym::FkResidualFuncCost2WrtFv1(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), encoder, p_a, p_b, p_c, p_d, sym::kDefaultEpsilon<double>);
                *H2 = (Matrix(4, 1) << Fkresidual_func_wrt_fv1).finished();
            }
            if(H3)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_Deltarot = sym::FkResidualFuncCost2WrtDeltarot(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), encoder, p_a, p_b, p_c, p_d, sym::kDefaultEpsilon<double>);
                *H3 = (Matrix(4, 3) << Fkresidual_func_wrt_Deltarot).finished();
            }
            if(H4)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_PositionVector = sym::FkResidualFuncCost2WrtPositionVector(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), encoder, p_a, p_b, p_c, p_d, sym::kDefaultEpsilon<double>);
                *H4 = (Matrix(4, 3) << Fkresidual_func_wrt_PositionVector).finished();
            }
            if(H5)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pa = sym::FkResidualFuncCost2WrtPa(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), encoder, p_a, p_b, p_c, p_d, sym::kDefaultEpsilon<double>);
                *H5 = (Matrix(4, 3) << Fkresidual_func_wrt_pa).finished();
            }
            if(H6)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pb = sym::FkResidualFuncCost2WrtPb(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), encoder, p_a, p_b, p_c, p_d, sym::kDefaultEpsilon<double>);
                *H6 = (Matrix(4, 3) << Fkresidual_func_wrt_pb).finished();
            }
            if(H7)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pc = sym::FkResidualFuncCost2WrtPc(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), encoder, p_a, p_b, p_c, p_d, sym::kDefaultEpsilon<double>);
                *H7 = (Matrix(4, 3) << Fkresidual_func_wrt_pc).finished();
            }
            if(H8)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pd = sym::FkResidualFuncCost2WrtPd(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), encoder, p_a, p_b, p_c, p_d, sym::kDefaultEpsilon<double>);
                *H8 = (Matrix(4, 3) << Fkresidual_func_wrt_pd).finished();
            }
            return (Vector(4) << Fkresidual_func).finished();
        }
    };


    class FK_factor_graoh_cost3 : public NoiseModelFactorN<double, double, gtsam::Rot3, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Point3>
    {

    private:
        Eigen::Vector4d encoder;    
        Eigen::Matrix3d rot_init;

    public:
        // Constructor
        FK_factor_graoh_cost3(Key key1, Key key2, Key key3, Key key4, Key key5, Key key6, Key key7, Key key8, Eigen::Vector4d encoder_, Eigen::Matrix3d rot_init_, const SharedNoiseModel &model) 
        : NoiseModelFactorN<double, double, gtsam::Rot3, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Point3>(model, key1, key2, key3, key4, key5, key6, key7, key8), encoder(encoder_), rot_init(rot_init_) {}

        // Evaluate the error
        Vector evaluateError(const double &fh1, const double &fv1, const gtsam::Rot3 &DeltaRot, const gtsam::Point3 &position_vector, const gtsam::Point3 &p_a, const gtsam::Point3 &p_b, const gtsam::Point3 &p_c, const gtsam::Point3 &p_d,
                             OptionalMatrixType H1,
                             OptionalMatrixType H2,
                             OptionalMatrixType H3,
                             OptionalMatrixType H4,
                             OptionalMatrixType H5,
                             OptionalMatrixType H6,
                             OptionalMatrixType H7,
                             OptionalMatrixType H8) const override
        {   
            Eigen::Matrix<double, 4, 1> Fkresidual_func = sym::FkResidualFuncCost3(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), encoder, p_a, p_b, p_c, p_d, sym::kDefaultEpsilon<double>);
            if(H1)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fh1 = sym::FkResidualFuncCost3WrtFh1(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), encoder, p_a, p_b, p_c, p_d, sym::kDefaultEpsilon<double>);
                *H1 = (Matrix(4, 1) << Fkresidual_func_wrt_fh1).finished();
            }
            if(H2)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fv1 = sym::FkResidualFuncCost3WrtFv1(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), encoder, p_a, p_b, p_c, p_d, sym::kDefaultEpsilon<double>);
                *H2 = (Matrix(4, 1) << Fkresidual_func_wrt_fv1).finished();
            }
            if(H3)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_Deltarot = sym::FkResidualFuncCost3WrtDeltarot(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), encoder, p_a, p_b, p_c, p_d, sym::kDefaultEpsilon<double>);
                *H3 = (Matrix(4, 3) << Fkresidual_func_wrt_Deltarot).finished();
            }
            if(H4)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_PositionVector = sym::FkResidualFuncCost3WrtPositionVector(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), encoder, p_a, p_b, p_c, p_d, sym::kDefaultEpsilon<double>);
                *H4 = (Matrix(4, 3) << Fkresidual_func_wrt_PositionVector).finished();
            }
            if(H5)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pa = sym::FkResidualFuncCost3WrtPa(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), encoder, p_a, p_b, p_c, p_d, sym::kDefaultEpsilon<double>);
                *H5 = (Matrix(4, 3) << Fkresidual_func_wrt_pa).finished();
            }
            if(H6)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pb = sym::FkResidualFuncCost3WrtPb(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), encoder, p_a, p_b, p_c, p_d, sym::kDefaultEpsilon<double>);
                *H6 = (Matrix(4, 3) << Fkresidual_func_wrt_pb).finished();
            }
            if(H7)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pc = sym::FkResidualFuncCost3WrtPc(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), encoder, p_a, p_b, p_c, p_d, sym::kDefaultEpsilon<double>);
                *H7 = (Matrix(4, 3) << Fkresidual_func_wrt_pc).finished();
            }
            if(H8)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pd = sym::FkResidualFuncCost3WrtPd(fh1, fv1, SymforceFromGtsam(DeltaRot), position_vector, SymforceFromGtsam(EigenMatrixToGtsamRot3(rot_init)), encoder, p_a, p_b, p_c, p_d, sym::kDefaultEpsilon<double>);
                *H8 = (Matrix(4, 3) << Fkresidual_func_wrt_pd).finished();
            }
            return (Vector(4) << Fkresidual_func).finished();
        }
    };

}