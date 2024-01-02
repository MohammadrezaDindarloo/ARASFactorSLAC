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
    class BetweenTwoPoint3FactorOneVariable : public NoiseModelFactor2<gtsam::Point3, double>
    {
        private:
        double measurement;
        gtsam::Point3 point2;
        
        public:
        // Constructor
        BetweenTwoPoint3FactorOneVariable(Key key1, Key key2, gtsam::Point3 point2_, double measurement_, const SharedNoiseModel &model) 
        : NoiseModelFactor2<gtsam::Point3, double>(model, key1, key2), point2(point2_), measurement(measurement_){}

        // Evaluate the error 
        Vector evaluateError(const gtsam::Point3 &point1, const double &offset,
                             OptionalMatrixType H1,
                             OptionalMatrixType H2) const override
        {
            if (H1)
            {   
                Eigen::Matrix<double, 1, 3> error_model_wrt_point1 = ErrorModelBetweenTowPointWrtPoint1(point1, point2, offset, measurement, sym::kDefaultEpsilon<double>);
                *H1 = (Matrix(1, 3) << error_model_wrt_point1).finished();
            }
            if (H2)
            {   
                Eigen::Matrix<double, 1, 1> error_model_wrt_ofset = ErrorModelBetweenTowPointWrtOfset(point1, point2, offset, measurement, sym::kDefaultEpsilon<double>);
                *H2 = (Matrix(1, 1) << error_model_wrt_ofset).finished();
            }
    
            Eigen::Matrix<double, 1, 1> error = sym::ErrorModelBetweenTowPoint(point1, point2, offset, measurement, sym::kDefaultEpsilon<double>) ;
            return (Vector(1) << error).finished();
        }
    };
}