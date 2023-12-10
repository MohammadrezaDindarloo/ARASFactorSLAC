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
    class FK_factor_graph_cost1 : public NoiseModelFactorN<double, double, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Point3>
    {

    private:     
        double rot_init_x; 
        double rot_init_y; 
        double rot_init_z; 
        double rot_init_w; 
        gtsam::Rot3 DeltaRot;
        gtsam::Point3 position_vector;

    public: 
        // Constructor
        FK_factor_graph_cost1(Key key1, Key key2, Key key3, Key key4, Key key5, Key key6, gtsam::Point3 position_vector_, gtsam::Rot3 DeltaRot_, double rot_init_x_, double rot_init_y_, double rot_init_z_, double rot_init_w_, const SharedNoiseModel &model) 
        : NoiseModelFactorN<double, double, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Point3>(model, key1, key2, key3, key4, key5, key6), position_vector(position_vector_), DeltaRot(DeltaRot_), rot_init_x(rot_init_x_), rot_init_y(rot_init_y_), rot_init_z(rot_init_z_), rot_init_w(rot_init_w_) {}

        // Evaluate the error
        Vector evaluateError(const double &fh1, const double &fv1, const gtsam::Point3 &p_a, const gtsam::Point3 &p_b, const gtsam::Point3 &p_c, const gtsam::Point3 &p_d,
                             OptionalMatrixType H1,
                             OptionalMatrixType H2,
                             OptionalMatrixType H3,
                             OptionalMatrixType H4,
                             OptionalMatrixType H5,
                             OptionalMatrixType H6) const override 
        {   
            Eigen::Matrix<double, 4, 1> Fkresidual_func = sym::FkResidualFuncCost1(fh1, fv1, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, sym::kDefaultEpsilon<double>);
            if(H1)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fh1 = sym::FkResidualFuncCost1WrtFh1(fh1, fv1, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, sym::kDefaultEpsilon<double>);
                *H1 = (Matrix(4, 1) << Fkresidual_func_wrt_fh1).finished();
            }
            if(H2)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fv1 = sym::FkResidualFuncCost1WrtFv1(fh1, fv1, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, sym::kDefaultEpsilon<double>);
                *H2 = (Matrix(4, 1) << Fkresidual_func_wrt_fv1).finished();
            }
            if(H3)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pa = sym::FkResidualFuncCost1WrtPa(fh1, fv1, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, sym::kDefaultEpsilon<double>);
                *H3 = (Matrix(4, 3) << Fkresidual_func_wrt_pa).finished();
            }
            if(H4)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pb = sym::FkResidualFuncCost1WrtPb(fh1, fv1, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, sym::kDefaultEpsilon<double>);
                *H4 = (Matrix(4, 3) << Fkresidual_func_wrt_pb).finished();
            }
            if(H5)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pc = sym::FkResidualFuncCost1WrtPc(fh1, fv1, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, sym::kDefaultEpsilon<double>);
                *H5 = (Matrix(4, 3) << Fkresidual_func_wrt_pc).finished();
            }
            if(H6)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pd = sym::FkResidualFuncCost1WrtPd(fh1, fv1, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, sym::kDefaultEpsilon<double>);
                *H6 = (Matrix(4, 3) << Fkresidual_func_wrt_pd).finished();
            }

            return (Vector(4) << Fkresidual_func).finished();
        }
    };


    class FK_factor_graph_cost2 : public NoiseModelFactorN<double, double, double, double, double, double, double, double, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Point3>
    {

    private:       
        double rot_init_x; 
        double rot_init_y; 
        double rot_init_z; 
        double rot_init_w; 
        gtsam::Rot3 DeltaRot;
        gtsam::Point3 position_vector;

    public:
        // Constructor
        FK_factor_graph_cost2(Key key1, Key key2, Key key3, Key key4, Key key5, Key key6, Key key7, Key key8, Key key9, Key key10, Key key11, Key key12, gtsam::Point3 position_vector_, gtsam::Rot3 DeltaRot_, double rot_init_x_, double rot_init_y_, double rot_init_z_, double rot_init_w_, const SharedNoiseModel &model) 
        : NoiseModelFactorN<double, double, double, double, double, double, double, double, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Point3>(model, key1, key2, key3, key4, key5, key6, key7, key8, key9, key10, key11, key12), position_vector(position_vector_), DeltaRot(DeltaRot_), rot_init_x(rot_init_x_), rot_init_y(rot_init_y_), rot_init_z(rot_init_z_), rot_init_w(rot_init_w_) {}

        // Evaluate the error
        Vector evaluateError(const double &fh1, const double &fv1, const double &fh2, const double &fv2, const double &fh3, const double &fv3, const double &fh4, const double &fv4, 
                             const gtsam::Point3 &p_a, const gtsam::Point3 &p_b, const gtsam::Point3 &p_c, const gtsam::Point3 &p_d,
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
            Eigen::Matrix<double, 4, 1> Fkresidual_func = sym::FkResidualFuncCost2(fh1, fv1, fh2, fv2, fh3, fv3, fh4, fv4, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, sym::kDefaultEpsilon<double>);
            if(H1)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fh1 = sym::FkResidualFuncCost2WrtFh1(fh1, fv1, fh2, fv2, fh3, fv3, fh4, fv4, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, sym::kDefaultEpsilon<double>);
                *H1 = (Matrix(4, 1) << Fkresidual_func_wrt_fh1).finished();
            }
            if(H2)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fv1 = sym::FkResidualFuncCost2WrtFv1(fh1, fv1, fh2, fv2, fh3, fv3, fh4, fv4, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, sym::kDefaultEpsilon<double>);
                *H2 = (Matrix(4, 1) << Fkresidual_func_wrt_fv1).finished();
            }
            if(H3)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fh2 = sym::FkResidualFuncCost2WrtFh2(fh1, fv1, fh2, fv2, fh3, fv3, fh4, fv4, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, sym::kDefaultEpsilon<double>);
                *H3 = (Matrix(4, 1) << Fkresidual_func_wrt_fh2).finished();
            }
            if(H4)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fv2 = sym::FkResidualFuncCost2WrtFv2(fh1, fv1, fh2, fv2, fh3, fv3, fh4, fv4, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, sym::kDefaultEpsilon<double>);
                *H4 = (Matrix(4, 1) << Fkresidual_func_wrt_fv2).finished();
            }
            if(H5)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fh3 = sym::FkResidualFuncCost2WrtFh3(fh1, fv1, fh2, fv2, fh3, fv3, fh4, fv4, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, sym::kDefaultEpsilon<double>);
                *H5 = (Matrix(4, 1) << Fkresidual_func_wrt_fh3).finished();
            }
            if(H6)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fv3 = sym::FkResidualFuncCost2WrtFv3(fh1, fv1, fh2, fv2, fh3, fv3, fh4, fv4, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, sym::kDefaultEpsilon<double>);
                *H6 = (Matrix(4, 1) << Fkresidual_func_wrt_fv3).finished();
            }
            if(H7)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fh4 = sym::FkResidualFuncCost2WrtFh4(fh1, fv1, fh2, fv2, fh3, fv3, fh4, fv4, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, sym::kDefaultEpsilon<double>);
                *H7 = (Matrix(4, 1) << Fkresidual_func_wrt_fh4).finished();
            }
            if(H8)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fv4 = sym::FkResidualFuncCost2WrtFv4(fh1, fv1, fh2, fv2, fh3, fv3, fh4, fv4, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, sym::kDefaultEpsilon<double>);
                *H8 = (Matrix(4, 1) << Fkresidual_func_wrt_fv4).finished();
            }
            if(H9)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pa = sym::FkResidualFuncCost2WrtPa(fh1, fv1, fh2, fv2, fh3, fv3, fh4, fv4, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, sym::kDefaultEpsilon<double>);
                *H9 = (Matrix(4, 3) << Fkresidual_func_wrt_pa).finished();
            }
            if(H10)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pb = sym::FkResidualFuncCost2WrtPb(fh1, fv1, fh2, fv2, fh3, fv3, fh4, fv4, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, sym::kDefaultEpsilon<double>);
                *H10 = (Matrix(4, 3) << Fkresidual_func_wrt_pb).finished();
            }
            if(H11)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pc = sym::FkResidualFuncCost2WrtPc(fh1, fv1, fh2, fv2, fh3, fv3, fh4, fv4, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, sym::kDefaultEpsilon<double>);
                *H11= (Matrix(4, 3) << Fkresidual_func_wrt_pc).finished();
            }
            if(H12)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pd = sym::FkResidualFuncCost2WrtPd(fh1, fv1, fh2, fv2, fh3, fv3, fh4, fv4, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, sym::kDefaultEpsilon<double>);
                *H12 = (Matrix(4, 3) << Fkresidual_func_wrt_pd).finished();
            }

            return (Vector(4) << Fkresidual_func).finished();
        }
    };


    class FK_factor_graph_cost3 : public NoiseModelFactorN<double, double, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Point3>
    {

    private:     
        double rot_init_x; 
        double rot_init_y; 
        double rot_init_z; 
        double rot_init_w; 
        gtsam::Rot3 DeltaRot;
        gtsam::Point3 position_vector;

    public:
        // Constructor
        FK_factor_graph_cost3(Key key1, Key key2, Key key3, Key key4, Key key5, Key key6, gtsam::Point3 position_vector_, gtsam::Rot3 DeltaRot_, double rot_init_x_, double rot_init_y_, double rot_init_z_, double rot_init_w_, const SharedNoiseModel &model) 
        : NoiseModelFactorN<double, double, gtsam::Point3, gtsam::Point3, gtsam::Point3, gtsam::Point3>(model, key1, key2, key3, key4, key5, key6), position_vector(position_vector_), DeltaRot(DeltaRot_), rot_init_x(rot_init_x_), rot_init_y(rot_init_y_), rot_init_z(rot_init_z_), rot_init_w(rot_init_w_) {}

        // Evaluate the error
        Vector evaluateError(const double &fh1, const double &fv1, const gtsam::Point3 &p_a, const gtsam::Point3 &p_b, const gtsam::Point3 &p_c, const gtsam::Point3 &p_d,
                             OptionalMatrixType H1,
                             OptionalMatrixType H2,
                             OptionalMatrixType H3,
                             OptionalMatrixType H4,
                             OptionalMatrixType H5,
                             OptionalMatrixType H6) const override       
        {   
            Eigen::Matrix<double, 4, 1> Fkresidual_func = sym::FkResidualFuncCost3(fh1, fv1, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, sym::kDefaultEpsilon<double>);
            if(H1)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fh1 = sym::FkResidualFuncCost3WrtFh1(fh1, fv1, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, sym::kDefaultEpsilon<double>);
                *H1 = (Matrix(4, 1) << Fkresidual_func_wrt_fh1).finished();
            }
            if(H2)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fv1 = sym::FkResidualFuncCost3WrtFv1(fh1, fv1, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, sym::kDefaultEpsilon<double>);
                *H2 = (Matrix(4, 1) << Fkresidual_func_wrt_fv1).finished();
            }
            if(H3)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pa = sym::FkResidualFuncCost3WrtPa(fh1, fv1, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, sym::kDefaultEpsilon<double>);
                *H3 = (Matrix(4, 3) << Fkresidual_func_wrt_pa).finished();
            }
            if(H4)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pb = sym::FkResidualFuncCost3WrtPb(fh1, fv1, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, sym::kDefaultEpsilon<double>);
                *H4 = (Matrix(4, 3) << Fkresidual_func_wrt_pb).finished();
            }
            if(H5)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pc = sym::FkResidualFuncCost3WrtPc(fh1, fv1, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, sym::kDefaultEpsilon<double>);
                *H5 = (Matrix(4, 3) << Fkresidual_func_wrt_pc).finished();
            }
            if(H6)
            {
                Eigen::Matrix<double, 4, 3> Fkresidual_func_wrt_pd = sym::FkResidualFuncCost3WrtPd(fh1, fv1, p_a, p_b, p_c, p_d, SymforceFromGtsam(DeltaRot), position_vector, rot_init_x, rot_init_y, rot_init_z, rot_init_w, sym::kDefaultEpsilon<double>);
                *H6 = (Matrix(4, 3) << Fkresidual_func_wrt_pd).finished();
            }

            return (Vector(4) << Fkresidual_func).finished();
        }
    };


    class FK_factor_graph_cost4 : public NoiseModelFactorN<double, double, double, double, double, double, double, double>
    {

    private:
        gtsam::Vector4 F_GT;

    public: 
        // Constructor
        FK_factor_graph_cost4(Key key1, Key key2, Key key3, Key key4, Key key5, Key key6, Key key7, Key key8, gtsam::Vector4 fc_, const SharedNoiseModel &model) 
        : NoiseModelFactorN<double, double, double, double, double, double, double, double>(model, key1, key2, key3, key4, key5, key6, key7, key8), F_GT(fc_) {}

        // Evaluate the error
        Vector evaluateError(const double &fh1, const double &fv1, const double &fh2, const double &fv2, const double &fh3, const double &fv3, const double &fh4, const double &fv4,
                             OptionalMatrixType H1,
                             OptionalMatrixType H2,
                             OptionalMatrixType H3,
                             OptionalMatrixType H4,
                             OptionalMatrixType H5,
                             OptionalMatrixType H6,
                             OptionalMatrixType H7,
                             OptionalMatrixType H8) const override 


        {   
            Eigen::Matrix<double, 4, 1> Fkresidual_func = sym::FkResidualFuncCost4(fh1, fv1, fh2, fv2, fh3, fv3, fh4, fv4, F_GT[0], F_GT[1], F_GT[2], F_GT[3], sym::kDefaultEpsilon<double>);
            if(H1)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fh1 = sym::FkResidualFuncCost4WrtFh1(fh1, fv1, fh2, fv2, fh3, fv3, fh4, fv4, F_GT[0], F_GT[1], F_GT[2], F_GT[3], sym::kDefaultEpsilon<double>);
                *H1 = (Matrix(4, 1) << Fkresidual_func_wrt_fh1).finished();
            }
            if(H2)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fv1 = sym::FkResidualFuncCost4WrtFv1(fh1, fv1, fh2, fv2, fh3, fv3, fh4, fv4, F_GT[0], F_GT[1], F_GT[2], F_GT[3], sym::kDefaultEpsilon<double>);
                *H2 = (Matrix(4, 1) << Fkresidual_func_wrt_fv1).finished();
            }
            if(H3)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fh2 = sym::FkResidualFuncCost4WrtFh2(fh1, fv1, fh2, fv2, fh3, fv3, fh4, fv4, F_GT[0], F_GT[1], F_GT[2], F_GT[3], sym::kDefaultEpsilon<double>);
                *H3 = (Matrix(4, 1) << Fkresidual_func_wrt_fh2).finished();
            }
            if(H4)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fv2 = sym::FkResidualFuncCost4WrtFv1(fh1, fv1, fh2, fv2, fh3, fv3, fh4, fv4, F_GT[0], F_GT[1], F_GT[2], F_GT[3], sym::kDefaultEpsilon<double>);
                *H4 = (Matrix(4, 1) << Fkresidual_func_wrt_fv2).finished();
            }
            if(H5)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fh3 = sym::FkResidualFuncCost4WrtFh3(fh1, fv1, fh2, fv2, fh3, fv3, fh4, fv4, F_GT[0], F_GT[1], F_GT[2], F_GT[3], sym::kDefaultEpsilon<double>);
                *H5 = (Matrix(4, 1) << Fkresidual_func_wrt_fh3).finished();
            }
            if(H6)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fv3 = sym::FkResidualFuncCost4WrtFv3(fh1, fv1, fh2, fv2, fh3, fv3, fh4, fv4, F_GT[0], F_GT[1], F_GT[2], F_GT[3], sym::kDefaultEpsilon<double>);
                *H6 = (Matrix(4, 1) << Fkresidual_func_wrt_fv3).finished();
            }
            if(H7)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fh4 = sym::FkResidualFuncCost4WrtFh4(fh1, fv1, fh2, fv2, fh3, fv3, fh4, fv4, F_GT[0], F_GT[1], F_GT[2], F_GT[3], sym::kDefaultEpsilon<double>);
                *H7 = (Matrix(4, 1) << Fkresidual_func_wrt_fh4).finished();
            }
            if(H8)
            {
                Eigen::Matrix<double, 4, 1> Fkresidual_func_wrt_fv4 = sym::FkResidualFuncCost4WrtFv4(fh1, fv1, fh2, fv2, fh3, fv3, fh4, fv4, F_GT[0], F_GT[1], F_GT[2], F_GT[3], sym::kDefaultEpsilon<double>);
                *H8 = (Matrix(4, 1) << Fkresidual_func_wrt_fv4).finished();
            }

            return (Vector(4) << Fkresidual_func).finished();
        }
    };


}
