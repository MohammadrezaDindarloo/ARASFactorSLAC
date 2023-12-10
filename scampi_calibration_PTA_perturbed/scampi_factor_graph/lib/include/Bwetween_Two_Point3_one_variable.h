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
    class BetweenTwoPoint3FactorOneVariable : public NoiseModelFactor1<gtsam::Point3>
    {
        private:
        Eigen::Matrix<double, 1, 1> measurement;
        gtsam::Point3 point2;
        
        public:
        // Constructor
        BetweenTwoPoint3FactorOneVariable(Key key1, gtsam::Point3 point2_, double measurement_, const SharedNoiseModel &model) 
        : NoiseModelFactor1<gtsam::Point3>(model, key1), point2(point2_), measurement(measurement_){}

        // Evaluate the error 
        Vector evaluateError(const gtsam::Point3 &point1,
                             OptionalMatrixType H1) const override
        {
            if (H1)
            {   
                Eigen::Matrix<double, 3, 1> error_model_wrt_point1 = ErrorModelWrtPoint1(point1, point2, sym::kDefaultEpsilon<double>);
                Eigen::Matrix<double, 1, 3> error_model_wrt_point1_ = error_model_wrt_point1.transpose();
                *H1 = (Matrix(1, 3) << error_model_wrt_point1_).finished();
            }
            
            Eigen::Matrix<double, 1, 1> error = sym::ErrorModel(point1, point2, sym::kDefaultEpsilon<double>) - measurement;
            return (Vector(1) << error).finished();
        }
    };
}