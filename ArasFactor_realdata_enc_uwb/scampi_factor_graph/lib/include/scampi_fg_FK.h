#pragma once

#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <boost/optional.hpp>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/inference/Symbol.h>
#include "scampi_function_header_include.h"
#include "scampi_fg_IK.h"


using namespace gtsam;
using namespace std;

namespace gtsam
{
    class FK_factor_graph_cost1 : public NoiseModelFactorN<double, double, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Pose3, gtsam::Rot3>
    {

    private:     

    public: 
        // Constructor
        FK_factor_graph_cost1(Key key1, Key key2, Key key3, Key key4, Key key5, Key key6, Key key7, Key key8, const SharedNoiseModel &model) 
        : NoiseModelFactorN<double, double, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Pose3, gtsam::Rot3>(model, key1, key2, key3, key4, key5, key6, key7, key8) {}

        // Evaluate the error
        Vector evaluateError(const double &fh1, const double &fv1, const gtsam::Point3 &p_a, const gtsam::Point3 &p_b, const gtsam::Point3 &p_c, const gtsam::Point3 &p_d, const  gtsam::Pose3 &TransformationMatrix, const  gtsam::Rot3 &DeltaRot,
                             OptionalMatrixType H1,
                             OptionalMatrixType H2,
                             OptionalMatrixType H3,
                             OptionalMatrixType H4,
                             OptionalMatrixType H5,
                             OptionalMatrixType H6,
                             OptionalMatrixType H7, 
                             OptionalMatrixType H8) const override 
        {  
            Eigen::Matrix<double, 4, 1> Fkresidual_func = sym::FkResidualFuncCost1(fh1, fv1, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
            if(H1)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fh1 = sym::FkResidualFuncCost1WrtFh1(fh1, fv1, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H1 = (Matrix(4, 1) << Fkresidual_func_wrt_fh1).finished();
            }
            if(H2)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fv1 = sym::FkResidualFuncCost1WrtFv1(fh1, fv1, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H2 = (Matrix(4, 1) << Fkresidual_func_wrt_fv1).finished();
            }
            if(H3)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pa = sym::FkResidualFuncCost1WrtPa(fh1, fv1, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H3 = (Matrix(4, 3) << Fkresidual_func_wrt_pa).finished();
            }
            if(H4)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pb = sym::FkResidualFuncCost1WrtPb(fh1, fv1, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H4 = (Matrix(4, 3) << Fkresidual_func_wrt_pb).finished();
            }
            if(H5)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pc = sym::FkResidualFuncCost1WrtPc(fh1, fv1, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H5 = (Matrix(4, 3) << Fkresidual_func_wrt_pc).finished();
            }
            if(H6)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pd = sym::FkResidualFuncCost1WrtPd(fh1, fv1, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H6 = (Matrix(4, 3) << Fkresidual_func_wrt_pd).finished();
            }
            if(H7)
            {
                Eigen::Matrix<double, 4, 6> Fkresidual_func_wrt_TransformationMatrix = sym::FkResidualFuncCost1WrtPose(fh1, fv1, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H7 = (Matrix(4, 6) << Fkresidual_func_wrt_TransformationMatrix).finished();
            }
            if(H8)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_Deltarot = sym::FkResidualFuncCost1WrtDeltarot(fh1, fv1, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H8 = (Matrix(4, 3) << Fkresidual_func_wrt_Deltarot).finished();
            }
            return (Vector(4) << Fkresidual_func).finished();
        }
    };


    class FK_factor_graph_cost2 : public NoiseModelFactorN<double, double, double, double, double, double, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Pose3, gtsam::Rot3>
    {

    private:       
        gtsam::Vector4 distance_measure;

    public:
        // Constructor
        FK_factor_graph_cost2(Key key1, Key key2, Key key3, Key key4, Key key5, Key key6, Key key7, Key key8, Key key9, Key key10, Key key11, Key key12, gtsam::Vector4 distance_measure_, const SharedNoiseModel &model) 
        : NoiseModelFactorN<double, double, double, double, double, double, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Pose3, gtsam::Rot3>(model, key1, key2, key3, key4, key5, key6, key7, key8, key9, key10, key11, key12), distance_measure(distance_measure_) {}

        // Evaluate the error
        Vector evaluateError(const double &fh1, const double &fv1, const double &ofset0, const double &ofset1, const double &ofset2, const double &ofset3,
                             const gtsam::Point3 &p_a, const gtsam::Point3 &p_b, const gtsam::Point3 &p_c, const gtsam::Point3 &p_d, const gtsam::Pose3 &TransformationMatrix, const gtsam::Rot3 &DeltaRot,
                             OptionalMatrixType H1,
                             OptionalMatrixType H2,
                             OptionalMatrixType H3,
                             OptionalMatrixType H4,
                             OptionalMatrixType H5,
                             OptionalMatrixType H6,
                             OptionalMatrixType H7,
                             OptionalMatrixType H8,
                             OptionalMatrixType H9,
                             OptionalMatrixType H10,
                             OptionalMatrixType H11,
                             OptionalMatrixType H12) const override
     
        { 
            Eigen::Matrix<double, 4, 1> Fkresidual_func = sym::FkResidualFuncCost2(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
            if(H1)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fh1 = sym::FkResidualFuncCost2WrtFh1(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H1 = (Matrix(4, 1) << Fkresidual_func_wrt_fh1).finished();
            }
            if(H2)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fv1 = sym::FkResidualFuncCost2WrtFv1(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H2 = (Matrix(4, 1) << Fkresidual_func_wrt_fv1).finished();
            }
            if(H3)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_ofset0 = sym::FkResidualFuncCost2WrtOfset0(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H3 = (Matrix(4, 1) << Fkresidual_func_wrt_ofset0).finished();
            }
            if(H4)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_ofset1 = sym::FkResidualFuncCost2WrtOfset1(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H4 = (Matrix(4, 1) << Fkresidual_func_wrt_ofset1).finished();
            }
            if(H5)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_ofset2 = sym::FkResidualFuncCost2WrtOfset2(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H5 = (Matrix(4, 1) << Fkresidual_func_wrt_ofset2).finished();
            }
            if(H6)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_ofset3 = sym::FkResidualFuncCost2WrtOfset3(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H6= (Matrix(4, 1) << Fkresidual_func_wrt_ofset3).finished();
            }
            if(H7)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pa= sym::FkResidualFuncCost2WrtPa(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H7 = (Matrix(4, 3) << Fkresidual_func_wrt_pa).finished();
            }
            if(H8)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pb = sym::FkResidualFuncCost2WrtPb(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H8 = (Matrix(4, 3) << Fkresidual_func_wrt_pb).finished();
            }
            if(H9)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pc = sym::FkResidualFuncCost2WrtPc(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H9 = (Matrix(4, 3) << Fkresidual_func_wrt_pc).finished();
            }
            if(H10)
            {              
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pd = sym::FkResidualFuncCost2WrtPd(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H10= (Matrix(4, 3) << Fkresidual_func_wrt_pd).finished();
            }
            if(H11)
            {
                Eigen::Matrix<double, 4, 6> Fkresidual_func_wrt_TransformationMatrix = sym::FkResidualFuncCost2WrtPose(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H11 = (Matrix(4, 6) << Fkresidual_func_wrt_TransformationMatrix).finished();
            }
            if(H12)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_Deltarot = sym::FkResidualFuncCost2WrtDeltarot(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H12 = (Matrix(4, 3) << Fkresidual_func_wrt_Deltarot).finished();
            }
            return (Vector(4) << Fkresidual_func).finished();
        }
    };


    class FK_factor_graph_cost3 : public NoiseModelFactorN<double, double, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Pose3, gtsam::Rot3>
    {

    private:     

    public:
        // Constructor
        FK_factor_graph_cost3(Key key1, Key key2, Key key3, Key key4, Key key5, Key key6, Key key7, Key key8, const SharedNoiseModel &model) 
        : NoiseModelFactorN<double, double, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Pose3, gtsam::Rot3>(model, key1, key2, key3, key4, key5, key6, key7, key8) {}

        // Evaluate the error
        Vector evaluateError(const double &fh1, const double &fv1, const gtsam::Point3 &p_a, const gtsam::Point3 &p_b, const gtsam::Point3 &p_c, const gtsam::Point3 &p_d, const gtsam::Pose3 &TransformationMatrix, const gtsam::Rot3 &DeltaRot,
                             OptionalMatrixType H1,
                             OptionalMatrixType H2,
                             OptionalMatrixType H3,
                             OptionalMatrixType H4,
                             OptionalMatrixType H5,
                             OptionalMatrixType H6,
                             OptionalMatrixType H7,
                             OptionalMatrixType H8) const override       
        {   
            Eigen::Matrix<double, 4, 1> Fkresidual_func = sym::FkResidualFuncCost3(fh1, fv1, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
            if(H1)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fh1 = sym::FkResidualFuncCost3WrtFh1(fh1, fv1, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H1 = (Matrix(4, 1) << Fkresidual_func_wrt_fh1).finished();
            }
            if(H2)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fv1 = sym::FkResidualFuncCost3WrtFv1(fh1, fv1, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H2 = (Matrix(4, 1) << Fkresidual_func_wrt_fv1).finished();
            }
            if(H3)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pa = sym::FkResidualFuncCost3WrtPa(fh1, fv1, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H3 = (Matrix(4, 3) << Fkresidual_func_wrt_pa).finished();
            }
            if(H4)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pb = sym::FkResidualFuncCost3WrtPb(fh1, fv1, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H4 = (Matrix(4, 3) << Fkresidual_func_wrt_pb).finished();
            }
            if(H5)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pc = sym::FkResidualFuncCost3WrtPc(fh1, fv1, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H5 = (Matrix(4, 3) << Fkresidual_func_wrt_pc).finished();
            }
            if(H6)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pd = sym::FkResidualFuncCost3WrtPd(fh1, fv1, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H6 = (Matrix(4, 3) << Fkresidual_func_wrt_pd).finished();
            }
            if(H7)
            {
                Eigen::Matrix<double, 4, 6> Fkresidual_func_wrt_TransformationMatrix = sym::FkResidualFuncCost3WrtPose(fh1, fv1, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H7 = (Matrix(4, 6) << Fkresidual_func_wrt_TransformationMatrix).finished();
            }
            if(H8)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_Deltarot = sym::FkResidualFuncCost3WrtDeltarot(fh1, fv1, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H8 = (Matrix(4, 3) << Fkresidual_func_wrt_Deltarot).finished();
            }
            return (Vector(4) << Fkresidual_func).finished();
        }
    };


    class FK_factor_graph_cost4 : public NoiseModelFactorN<double, double, double, double, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Point3, double, double, gtsam::Pose3, gtsam::Rot3>
    {

    private:
        gtsam::Vector4 distance_measure;

    public: 
        // Constructor
        FK_factor_graph_cost4(Key key1, Key key2, Key key3, Key key4, Key key5, Key key6, Key key7, Key key8, Key key9, Key key10, Key key11, Key key12, gtsam::Vector4 distance_measure_, const SharedNoiseModel &model) 
        : NoiseModelFactorN<double, double, double, double, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Point3, double, double, gtsam::Pose3, gtsam::Rot3>(model, key1, key2, key3, key4, key5, key6, key7, key8, key9, key10, key11, key12), distance_measure(distance_measure_) {}

        // Evaluate the error
        Vector evaluateError(const double &ofset0, const double &ofset1, const double &ofset2, const double &ofset3, const gtsam::Point3 &p_a, const gtsam::Point3 &p_b, const gtsam::Point3 &p_c, const gtsam::Point3 &p_d, const double &fh1, const double &fv1, const gtsam::Pose3 &TransformationMatrix, const gtsam::Rot3 &DeltaRot,
                             OptionalMatrixType H1,
                             OptionalMatrixType H2,
                             OptionalMatrixType H3,
                             OptionalMatrixType H4,
                             OptionalMatrixType H5,
                             OptionalMatrixType H6,
                             OptionalMatrixType H7,
                             OptionalMatrixType H8,
                             OptionalMatrixType H9,
                             OptionalMatrixType H10,
                             OptionalMatrixType H11,
                             OptionalMatrixType H12) const override 

        {   
            Eigen::Matrix<double, 4, 1> Fkresidual_func = sym::FkResidualFuncCost4(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
            if(H1)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_ofset0 = sym::FkResidualFuncCost4WrtOfset0(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H1 = (Matrix(4, 1) << Fkresidual_func_wrt_ofset0).finished();
            }
            if(H2)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_ofset1 = sym::FkResidualFuncCost4WrtOfset1(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H2 = (Matrix(4, 1) << Fkresidual_func_wrt_ofset1).finished();
            }
            if(H3)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_ofset2 = sym::FkResidualFuncCost4WrtOfset2(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H3 = (Matrix(4, 1) << Fkresidual_func_wrt_ofset2).finished();
            }
            if(H4)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_ofset3 = sym::FkResidualFuncCost4WrtOfset3(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H4 = (Matrix(4, 1) << Fkresidual_func_wrt_ofset3).finished();
            }
            if(H5)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pa = sym::FkResidualFuncCost4WrtPa(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H5 = (Matrix(4, 3) << Fkresidual_func_wrt_pa).finished();
            }
            if(H6)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pb = sym::FkResidualFuncCost4WrtPb(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H6 = (Matrix(4, 3) << Fkresidual_func_wrt_pb).finished();
            }
            if(H7)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pc = sym::FkResidualFuncCost4WrtPc(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H7 = (Matrix(4, 3) << Fkresidual_func_wrt_pc).finished();
            }
            if(H8)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pd = sym::FkResidualFuncCost4WrtPd(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H8 = (Matrix(4, 3) << Fkresidual_func_wrt_pd).finished();
            }
            if(H9)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fh1 = sym::FkResidualFuncCost4WrtFh1(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H9 = (Matrix(4, 1) << Fkresidual_func_wrt_fh1).finished();
            }
            if(H10)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fv1 = sym::FkResidualFuncCost4WrtFv1(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H10 = (Matrix(4, 1) << Fkresidual_func_wrt_fv1).finished();
            }
            if(H11)
            {
                Eigen::Matrix<double, 4, 6> Fkresidual_func_wrt_TransformationMatrix = sym::FkResidualFuncCost4WrtPose(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H11 = (Matrix(4, 6) << Fkresidual_func_wrt_TransformationMatrix).finished();
            }
            if(H12)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_Deltarot = sym::FkResidualFuncCost4WrtDeltarot(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H12 = (Matrix(4, 3) << Fkresidual_func_wrt_Deltarot).finished();
            }
            return (Vector(4) << Fkresidual_func).finished();
        }
    };



    class FK_factor_graph_cost5 : public NoiseModelFactorN<double, double, double, double, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Point3, double, double, gtsam::Pose3, gtsam::Rot3>
    {

    private:
        gtsam::Vector4 distance_measure;

    public: 
        // Constructor
        FK_factor_graph_cost5(Key key1, Key key2, Key key3, Key key4, Key key5, Key key6, Key key7, Key key8, Key key9, Key key10, Key key11, Key key12, gtsam::Vector4 distance_measure_, const SharedNoiseModel &model) 
        : NoiseModelFactorN<double, double, double, double, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Point3, double, double, gtsam::Pose3, gtsam::Rot3>(model, key1, key2, key3, key4, key5, key6, key7, key8, key9, key10, key11, key12), distance_measure(distance_measure_) {}

        // Evaluate the error
        Vector evaluateError(const double &ofset0, const double &ofset1, const double &ofset2, const double &ofset3, const gtsam::Point3 &p_a, const gtsam::Point3 &p_b, const gtsam::Point3 &p_c, const gtsam::Point3 &p_d, const double &fh1, const double &fv1, const gtsam::Pose3 &TransformationMatrix, const gtsam::Rot3 &DeltaRot,
                             OptionalMatrixType H1,
                             OptionalMatrixType H2,
                             OptionalMatrixType H3,
                             OptionalMatrixType H4,
                             OptionalMatrixType H5,
                             OptionalMatrixType H6,
                             OptionalMatrixType H7,
                             OptionalMatrixType H8,
                             OptionalMatrixType H9,
                             OptionalMatrixType H10,
                             OptionalMatrixType H11,
                             OptionalMatrixType H12) const override 

        {   
            Eigen::Matrix<double, 4, 1> Fkresidual_func = sym::FkResidualFuncCost5(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
            if(H1)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_ofset0 = sym::FkResidualFuncCost5WrtOfset0(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H1 = (Matrix(4, 1) << Fkresidual_func_wrt_ofset0).finished();
            }
            if(H2)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_ofset1 = sym::FkResidualFuncCost5WrtOfset1(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H2 = (Matrix(4, 1) << Fkresidual_func_wrt_ofset1).finished();
            }
            if(H3)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_ofset2 = sym::FkResidualFuncCost5WrtOfset2(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H3 = (Matrix(4, 1) << Fkresidual_func_wrt_ofset2).finished();
            }
            if(H4)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_ofset3 = sym::FkResidualFuncCost5WrtOfset3(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H4 = (Matrix(4, 1) << Fkresidual_func_wrt_ofset3).finished();
            }
            if(H5)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pa = sym::FkResidualFuncCost5WrtPa(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H5 = (Matrix(4, 3) << Fkresidual_func_wrt_pa).finished();
            }
            if(H6)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pb = sym::FkResidualFuncCost5WrtPb(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H6 = (Matrix(4, 3) << Fkresidual_func_wrt_pb).finished();
            }
            if(H7)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pc = sym::FkResidualFuncCost5WrtPc(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H7 = (Matrix(4, 3) << Fkresidual_func_wrt_pc).finished();
            }
            if(H8)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pd = sym::FkResidualFuncCost5WrtPd(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H8 = (Matrix(4, 3) << Fkresidual_func_wrt_pd).finished();
            }
            if(H9)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fh1 = sym::FkResidualFuncCost5WrtFh1(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H9 = (Matrix(4, 1) << Fkresidual_func_wrt_fh1).finished();
            }
            if(H10)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fv1 = sym::FkResidualFuncCost5WrtFv1(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H10 = (Matrix(4, 1) << Fkresidual_func_wrt_fv1).finished();
            }
            if(H11)
            {
                Eigen::Matrix<double, 4, 6> Fkresidual_func_wrt_TransformationMatrix = sym::FkResidualFuncCost5WrtPose(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H11 = (Matrix(4, 6) << Fkresidual_func_wrt_TransformationMatrix).finished();
            }
            if(H12)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_Deltarot = sym::FkResidualFuncCost5WrtDeltarot(fh1, fv1, ofset0, ofset1, ofset2, ofset3, Eigen::Matrix(distance_measure), p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), SymPoseFromGtsam(TransformationMatrix), sym::kDefaultEpsilon<double>);
                *H12 = (Matrix(4, 3) << Fkresidual_func_wrt_Deltarot).finished();
            }
            return (Vector(4) << Fkresidual_func).finished();
        }
    };
}