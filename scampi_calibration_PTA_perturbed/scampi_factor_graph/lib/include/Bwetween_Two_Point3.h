#pragma once

#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <boost/optional.hpp>
#include <gtsam/inference/Symbol.h>
#include <gtsam/geometry/Point3.h>
#include "scampi_function_header_include.h"



using namespace gtsam;
using namespace std;
using namespace sym;

namespace gtsam
{   
    class BetweenTwoPoint3Factor : public NoiseModelFactor2<gtsam::Point3, gtsam::Point3>
    {
        private:
        Eigen::Matrix<double, 1, 1> measurement;
        
        public:
        // Constructor
        BetweenTwoPoint3Factor(Key key1, Key key2, double measurement_, const SharedNoiseModel &model) 
        : NoiseModelFactor2<gtsam::Point3, gtsam::Point3>(model, key1, key2), measurement(measurement_){}

        // Evaluate the error 
        Vector evaluateError(const gtsam::Point3 &point1, const gtsam::Point3 &point2,
                             OptionalMatrixType H1, OptionalMatrixType H2) const override
        {
            if (H1)
            {   
                Eigen::Matrix<double, 3, 1> error_model_wrt_point1 = ErrorModelWrtPoint1(point1, point2, sym::kDefaultEpsilon<double>);
                Eigen::Matrix<double, 1, 3> error_model_wrt_point1_ = error_model_wrt_point1.transpose();
                *H1 = (Matrix(1, 3) << error_model_wrt_point1_).finished();
            }
            if (H2)
            {   
                Eigen::Matrix<double, 3, 1> error_model_wrt_point2 = ErrorModelWrtPoint2(point1, point2, sym::kDefaultEpsilon<double>);
                Eigen::Matrix<double, 1, 3> error_model_wrt_point2_ = error_model_wrt_point2.transpose();
                *H2 = (Matrix(1, 3) << error_model_wrt_point2_).finished();
            }
            
            Eigen::Matrix<double, 1, 1> error = sym::ErrorModel(point1, point2, sym::kDefaultEpsilon<double>) - measurement;
            return (Vector(1) << error).finished();
        }
    };
}