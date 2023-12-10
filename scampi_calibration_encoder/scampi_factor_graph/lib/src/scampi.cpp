#include "../include/libscampi.h"
#include <gtsam/slam/PriorFactor.h>

// ****************************************** IK optimization *******************************************
// a function to create ivnerse kinematic factor graph 
void inverse_kinematic_factor_graph_optimizer(double p_init_0, double p_init_1, double p_init_2,
                              double rot_init_x, double rot_init_y, double rot_init_z, double rot_init_w, 
                              RobotParameters<double> params,
                              int largest_cable,
                              double init_estimate_h1, double init_estimate_v1, gtsam::Rot3 init_estimate_rot, gtsam::Values *oprimization_result_LM)
{
    NonlinearFactorGraph graph;
    Values initial_estimate;

    auto Sensor_noiseModel_cost1 = gtsam::noiseModel::Isotropic::Sigma(4, sqrt(500));
    auto Sensor_noiseModel_cost2 = gtsam::noiseModel::Isotropic::Sigma(4, sqrt(10));
    auto Sensor_noiseModel_cost3 = gtsam::noiseModel::Isotropic::Sigma(4, 1);

    graph.add(std::make_shared<IK_factor_graoh_cost1>(Symbol('h', 1), Symbol('v', 1), Symbol('r', 1), params.pulleys[0], params.pulleys[1], params.pulleys[2], params.pulleys[3], p_init_0, p_init_1, p_init_2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, largest_cable, Sensor_noiseModel_cost1));
    graph.add(std::make_shared<IK_factor_graoh_cost2>(Symbol('h', 1), Symbol('v', 1), Symbol('r', 1), params.pulleys[0], params.pulleys[1], params.pulleys[2], params.pulleys[3], p_init_0, p_init_1, p_init_2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, largest_cable, Sensor_noiseModel_cost2));
    graph.add(std::make_shared<IK_factor_graoh_cost3>(Symbol('h', 1), Symbol('v', 1), Symbol('r', 1), params.pulleys[0], params.pulleys[1], params.pulleys[2], params.pulleys[3], p_init_0, p_init_1, p_init_2, rot_init_x, rot_init_y, rot_init_z, rot_init_w, largest_cable, Sensor_noiseModel_cost3));

    initial_estimate.insert(Symbol('h', 1), init_estimate_h1);
    initial_estimate.insert(Symbol('v', 1), init_estimate_v1);
    initial_estimate.insert(Symbol('r', 1), init_estimate_rot);

    LevenbergMarquardtOptimizer optimizer(graph, initial_estimate);
    Values result_LM = optimizer.optimize();
    std::cout << std::endl << "-------------------catenary result--------------------------" << std::endl;
    std::cout << std::endl << "inverse kinematic optimization error: " << optimizer.error() << std::endl;
    *oprimization_result_LM = result_LM;
}

// a function to hold parameters and invoke optimizer and back the optimized data
void ikSolver(RobotParameters<double> params, 
              Eigen::Matrix3d rot_init, 
              Eigen::Vector3d p_platform,
              IKDataOut<double> *result)
{
    //reorder the cable forces and choose the cable with largest length as the the first cable (for numerical stability)
    VectorXi reorder_idx(params.pulleys.size());
    RobotParameters<double> params_reord;
    RobotState<double> state;
    state.rot_platform = rot_init;
    state.p_platform = p_platform;
    changeOrderForSolver<double>(state, params, &params_reord, &reorder_idx);
    //Compute initil cable forces as starting points for the solver
    double fh0, fv0;
    computeInitCableForces<double>(&fh0, &fv0, p_platform, rot_init, params_reord);
    // Initialize the quaternion array
    double rotation_matrix_double[9];
    for (int i = 0; i < 9; ++i) {
        rotation_matrix_double[i] = rot_init.data()[i];
    }
    double quaternion[4];    
    ceres::RotationMatrixToQuaternion<double>(rotation_matrix_double, quaternion);
    double rot_init_x = quaternion[1];
    double rot_init_y = quaternion[2];
    double rot_init_z = quaternion[3];
    double rot_init_w = quaternion[0];
    // initial p_platfrom
    double p_init_0 = p_platform[0];
    double p_init_1 = p_platform[1];
    double p_init_2 = p_platform[2];
    // define largest cable
    int largest_cable = -1;
    if(reorder_idx[0] == 0 && reorder_idx[1] == 1 && reorder_idx[2] == 2 && reorder_idx[3] == 3)
    {
        largest_cable = 0; 
    }
    else if (reorder_idx[0] == 0 && reorder_idx[1] == 1 && reorder_idx[2] == 3 && reorder_idx[3] == 2)
    {
        largest_cable = 1; 
    }
    else if (reorder_idx[0] == 0 && reorder_idx[1] == 2 && reorder_idx[2] == 1 && reorder_idx[3] == 3)
    {
        largest_cable = 2; 
    }
    else if (reorder_idx[0] == 0 && reorder_idx[1] == 2 && reorder_idx[2] == 3 && reorder_idx[3] == 1)
    {
        largest_cable = 3; 
    }
    else if (reorder_idx[0] == 0 && reorder_idx[1] == 3 && reorder_idx[2] == 1 && reorder_idx[3] == 2)
    {
        largest_cable = 4; 
    }
    else if (reorder_idx[0] == 0 && reorder_idx[1] == 3 && reorder_idx[2] == 2 && reorder_idx[3] == 1)
    {
        largest_cable = 5; 
    }
    else if (reorder_idx[0] == 1 && reorder_idx[1] == 0 && reorder_idx[2] == 2 && reorder_idx[3] == 3)
    {
        largest_cable = 6; 
    }
    else if (reorder_idx[0] == 1 && reorder_idx[1] == 0 && reorder_idx[2] == 3 && reorder_idx[3] == 2)
    {
        largest_cable = 7; 
    }
    else if (reorder_idx[0] == 1 && reorder_idx[1] == 2 && reorder_idx[2] == 0 && reorder_idx[3] == 3)
    {
        largest_cable = 8; 
    }
    else if (reorder_idx[0] == 1 && reorder_idx[1] == 2 && reorder_idx[2] == 3 && reorder_idx[3] == 0)
    {
        largest_cable = 9; 
    }
    else if (reorder_idx[0] == 1 && reorder_idx[1] == 3 && reorder_idx[2] == 0 && reorder_idx[3] == 2)
    {
        largest_cable = 10; 
    }
    else if (reorder_idx[0] == 1 && reorder_idx[1] == 3 && reorder_idx[2] == 2 && reorder_idx[3] == 0)
    {
        largest_cable = 11; 
    }
    else if (reorder_idx[0] == 2 && reorder_idx[1] == 0 && reorder_idx[2] == 1 && reorder_idx[3] == 3)
    {
        largest_cable = 12; 
    }
    else if (reorder_idx[0] == 2 && reorder_idx[1] == 0 && reorder_idx[2] == 3 && reorder_idx[3] == 1)
    {
        largest_cable = 13; 
    }
    else if (reorder_idx[0] == 2 && reorder_idx[1] == 1 && reorder_idx[2] == 0 && reorder_idx[3] == 3)
    {
        largest_cable = 14; 
    }
    else if (reorder_idx[0] == 2 && reorder_idx[1] == 1 && reorder_idx[2] == 3 && reorder_idx[3] == 0)
    {
        largest_cable = 15; 
    }
    else if (reorder_idx[0] == 2 && reorder_idx[1] == 3 && reorder_idx[2] == 0 && reorder_idx[3] == 1)
    {
        largest_cable = 16; 
    }
    else if (reorder_idx[0] == 2 && reorder_idx[1] == 3 && reorder_idx[2] == 1 && reorder_idx[3] == 0)
    {
        largest_cable = 17; 
    }
    else if (reorder_idx[0] == 3 && reorder_idx[1] == 0 && reorder_idx[2] == 1 && reorder_idx[3] == 2)
    {
        largest_cable = 18; 
    }
    else if (reorder_idx[0] == 3 && reorder_idx[1] == 0 && reorder_idx[2] == 2 && reorder_idx[3] == 1)
    {
        largest_cable = 19; 
    }
    else if (reorder_idx[0] == 3 && reorder_idx[1] == 1 && reorder_idx[2] == 0 && reorder_idx[3] == 2)
    {
        largest_cable = 20; 
    }
    else if (reorder_idx[0] == 3 && reorder_idx[1] == 1 && reorder_idx[2] == 2 && reorder_idx[3] == 0)
    {
        largest_cable = 21; 
    }
    else if (reorder_idx[0] == 3 && reorder_idx[1] == 2 && reorder_idx[2] == 0 && reorder_idx[3] == 1)
    {
        largest_cable = 22; 
    }
    else if (reorder_idx[0] == 3 && reorder_idx[1] == 2 && reorder_idx[2] == 1 && reorder_idx[3] == 0)
    {
        largest_cable = 23; 
    }
    else
    {
        std::cout << "Error: Cable index is wrong!!" << std::endl;
        exit(1);
    }
    // initial values for variable 
    double init_estimate_h1 = fh0;
    double init_estimate_v1 = -fv0;
    gtsam::Rot3 init_estimate_rot = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    // run optimization!
    gtsam::Values optimization_result;
    inverse_kinematic_factor_graph_optimizer(p_init_0, p_init_1, p_init_2,
                                             rot_init_x, rot_init_y, rot_init_z, rot_init_w, 
                                             params,
                                             largest_cable,
                                             init_estimate_h1, init_estimate_v1, init_estimate_rot, &optimization_result);

    // optimization_result.print();
    // harvest the results
    double fh = optimization_result.at<double>(Symbol('h', 1)); //optimized horizontal force for the first cable
    double fv = optimization_result.at<double>(Symbol('v', 1)); //optimized vertical force for the first cable
    // Extract the optimized orientation matrix of the moving platform
    gtsam::Rot3 rot_optimized = optimization_result.at<gtsam::Rot3>(Symbol('r', 1));
    Eigen::Matrix3d rot_optimized_ = gtsamRot3ToEigenMatrix(rot_optimized);
    Matrix3d rot_result = rot_init * rot_optimized_;
    result->rot_platform = rot_result;

    // Use the utils functions once again to compute the force in other cables and the catenary variables
    GeometricVariables<double> geom_vars;
    CatenaryVariables<double> cat_vars;

    state.rot_platform = rot_result; // Now we know the optimzed value of end-effector position so just use it

    getGeometricVariables<double>(state,params_reord,&geom_vars);
    getCableForces<double>(fh, fv, &state, params_reord,geom_vars);
    getCatenaryVariables<double>(state,params_reord, geom_vars,&cat_vars);
    //reverse the order of cables back to the normal configuration (base on notebook index)
    reverseOrderForSolver<double>(state, geom_vars, cat_vars, result, reorder_idx);
}

// ****************************************** FK optimization *******************************************
// a function to create forward kinematic factor graph 
void forward_kinematic_factor_graph_optimizer(std::vector<Eigen::Matrix<double, 2, 1>> cable_forces_collection , RobotParameters<double> params, 
                              std::vector<std::vector< Eigen::Matrix<double, 3, 1> >> b_in_w_collection, std::vector<std::vector< Eigen::Matrix<double, 3, 1> >> p_in_w_collection,
                              std::vector<Eigen::Matrix<double, 3, 1>> pulley_position_estimate,
                              std::vector<Eigen::Matrix<double, 4, 1>> quaternion_rot_init_collection, 
                              gtsam::Rot3 DeltaRot, std::vector<Eigen::Matrix<double, 3, 1>> p_platform_collection, gtsam::Values *oprimization_result_LM)
{
    NonlinearFactorGraph graph;
    Values initial_estimate;

    auto Sensor_noiseModel_cost1 = gtsam::noiseModel::Isotropic::Sigma(4, sqrt(10));
    auto Sensor_noiseModel_cost2 = gtsam::noiseModel::Isotropic::Sigma(4, sqrt(10));
    auto Sensor_noiseModel_cost3 = gtsam::noiseModel::Isotropic::Sigma(4, sqrt(1));
    auto Sensor_noiseModel_cost4 = gtsam::noiseModel::Isotropic::Sigma(4, sqrt(10));
    noiseModel::Isotropic::shared_ptr Sensor_noiseModel_betweenfactor = gtsam::noiseModel::Isotropic::Sigma(1, sqrt(0.01));
    noiseModel::Isotropic::shared_ptr Sensor_noiseModel_betweenfactor_one_variable = gtsam::noiseModel::Isotropic::Sigma(1, sqrt(0.01));
    noiseModel::Isotropic::shared_ptr prior_noiseModel = gtsam::noiseModel::Isotropic::Sigma(1, sqrt(0.01));
    noiseModel::Isotropic::shared_ptr prior_noiseModel_3D = gtsam::noiseModel::Isotropic::Sigma(3, sqrt(0.01));

    for (size_t i = 0; i < p_platform_collection.size(); i++)
    {
        double b1_p1_distance = (b_in_w_collection[i][0] - p_in_w_collection[i][0]).norm()+1.1;
        double b2_p2_distance = (b_in_w_collection[i][1] - p_in_w_collection[i][1]).norm()-0.5;
        double b3_p3_distance = (b_in_w_collection[i][2] - p_in_w_collection[i][2]).norm()+2;
        double b4_p4_distance = (b_in_w_collection[i][3] - p_in_w_collection[i][3]).norm()+1.8;
        gtsam::Vector4 distance_measure = {b1_p1_distance, b2_p2_distance, b3_p3_distance, b4_p4_distance};

        // graph.add(std::make_shared<FK_factor_graph_cost1>(Symbol('h', i), Symbol('v', i), Symbol('p', 1), Symbol('p', 2), Symbol('p', 3), Symbol('p', 4), p_platform_collection[i], DeltaRot, quaternion_rot_init_collection[i][0], quaternion_rot_init_collection[i][1], quaternion_rot_init_collection[i][2], quaternion_rot_init_collection[i][3], Sensor_noiseModel_cost1));
        // graph.add(std::make_shared<FK_factor_graph_cost2>(Symbol('h', i), Symbol('v', i), Symbol('o', 0), Symbol('o', 1), Symbol('o', 2), Symbol('o', 3), Symbol('p', 1), Symbol('p', 2), Symbol('p', 3), Symbol('p', 4), distance_measure, p_platform_collection[i], DeltaRot, quaternion_rot_init_collection[i][0], quaternion_rot_init_collection[i][1], quaternion_rot_init_collection[i][2], quaternion_rot_init_collection[i][3], Sensor_noiseModel_cost2));
        // graph.add(std::make_shared<FK_factor_graph_cost3>(Symbol('h', i), Symbol('v', i), Symbol('p', 1), Symbol('p', 2), Symbol('p', 3), Symbol('p', 4), p_platform_collection[i], DeltaRot, quaternion_rot_init_collection[i][0], quaternion_rot_init_collection[i][1], quaternion_rot_init_collection[i][2], quaternion_rot_init_collection[i][3], Sensor_noiseModel_cost3));
        graph.add(std::make_shared<FK_factor_graph_cost4>(Symbol('o', 0), Symbol('o', 1), Symbol('o', 2), Symbol('o', 3), Symbol('p', 1), Symbol('p', 2), Symbol('p', 3), Symbol('p', 4), Symbol('h', i), Symbol('v', i), distance_measure, p_platform_collection[i], DeltaRot, quaternion_rot_init_collection[i][0], quaternion_rot_init_collection[i][1], quaternion_rot_init_collection[i][2], quaternion_rot_init_collection[i][3], Sensor_noiseModel_cost4));
    }

    // initial_estimate for variables
    initial_estimate.insert(Symbol('p', 1), pulley_position_estimate[0]);
    initial_estimate.insert(Symbol('p', 2), pulley_position_estimate[1]);
    initial_estimate.insert(Symbol('p', 3), pulley_position_estimate[2]);
    initial_estimate.insert(Symbol('p', 4), pulley_position_estimate[3]);

    initial_estimate.insert(Symbol('o', 0), 0.0);
    initial_estimate.insert(Symbol('o', 1), 0.0);
    initial_estimate.insert(Symbol('o', 2), 0.0);
    initial_estimate.insert(Symbol('o', 3), 0.0);
    for (size_t i = 0; i < p_platform_collection.size(); i++)
    {
        initial_estimate.insert(Symbol('h', i), cable_forces_collection[i].norm()/sqrt(2));
        initial_estimate.insert(Symbol('v', i), cable_forces_collection[i].norm()/sqrt(2));   
    }

    LevenbergMarquardtOptimizer optimizer(graph, initial_estimate);
    Values result_LM = optimizer.optimize();
    std::cout << std::endl << "-------------------catenary result--------------------------" << std::endl;
    std::cout << std::endl << "forward kinematic optimization error: " << optimizer.error() << std::endl;
    std::cout << std::endl << "ofset 0: " << result_LM.at<double>(Symbol('o', 0)) << std::endl;
    std::cout << std::endl << "ofset 1: " << result_LM.at<double>(Symbol('o', 1)) << std::endl;
    std::cout << std::endl << "ofset 2: " << result_LM.at<double>(Symbol('o', 2)) << std::endl;
    std::cout << std::endl << "ofset 3: " << result_LM.at<double>(Symbol('o', 3)) << std::endl;

    *oprimization_result_LM = result_LM;
}

// a function to hold parameters and invoke optimizer and back the optimized data
void fkSolver(std::vector<Eigen::Matrix<double, 3, 1>> p_platform_collection,  
              std::vector<Eigen::Matrix<double, 3, 3>> rot_platform_collection, 
              std::vector<Eigen::Matrix<double, 3, 1>> pulley_position_estimate,
              std::vector<Eigen::Matrix<double, 2, 1>> cable_forces_collection, 
              RobotParameters<double> params, 
              FKDataOut<double> *result)
{   
    RobotState<double> state_init;
    GeometricVariables<double> geom_vars_calibration;
    std::vector<Eigen::Matrix<double, 4, 1>> quaternion_rot_init_collection;
    std::vector<std::vector< Eigen::Matrix<double, 3, 1> >> b_in_w_collection;
    std::vector<std::vector< Eigen::Matrix<double, 3, 1> >> p_in_w_collection;
    for (size_t j = 0; j < p_platform_collection.size(); j++)
    {
        // Initialize the quaternion array
        double rotation_matrix_double[9];
        for (int i = 0; i < 9; ++i) {
            rotation_matrix_double[i] = rot_platform_collection[j].data()[i];
        }
        double quaternion[4];    
        ceres::RotationMatrixToQuaternion<double>(rotation_matrix_double, quaternion);
        quaternion_rot_init_collection.push_back(Eigen::Matrix<double, 4, 1>(quaternion[1], quaternion[2], quaternion[3], quaternion[0]));
        state_init.p_platform = p_platform_collection[j];
        state_init.rot_platform = rot_platform_collection[j];
        getGeometricVariables<double>(state_init,params,&geom_vars_calibration);
        b_in_w_collection.push_back(geom_vars_calibration.b_in_w);
        p_in_w_collection.push_back(geom_vars_calibration.p_in_w);
    }
    // initial values for variable 
    gtsam::Rot3 init_estimate_rotation = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    // run optimization!
    gtsam::Values optimization_result; 
    forward_kinematic_factor_graph_optimizer(cable_forces_collection, params, b_in_w_collection, p_in_w_collection,
                                             pulley_position_estimate, quaternion_rot_init_collection, 
                                             init_estimate_rotation, p_platform_collection, &optimization_result);

    // optimization_result.print();
    double disparity_in_f_h = 0.0;
    double disparity_in_f_v = 0.0;
    for (size_t i = 0; i < cable_forces_collection.size(); i++)
    {
        double temp1 = abs(optimization_result.at<double>(Symbol('h', i)) - cable_forces_collection[i][0]);
        double temp2 = abs(optimization_result.at<double>(Symbol('v', i)) - cable_forces_collection[i][1]);
        disparity_in_f_h = disparity_in_f_h + temp1;
        disparity_in_f_v = disparity_in_f_v + temp2;
    }
    std::cout << std::endl << "disparity_in_f_h: " << disparity_in_f_h << std::endl;
    std::cout << std::endl << "disparity_in_f_v: " << disparity_in_f_v << std::endl;
    //harvest the results
    result->pulley_position.push_back(optimization_result.at<gtsam::Point3>(Symbol('p', 1)));
    result->pulley_position.push_back(optimization_result.at<gtsam::Point3>(Symbol('p', 2)));
    result->pulley_position.push_back(optimization_result.at<gtsam::Point3>(Symbol('p', 3)));
    result->pulley_position.push_back(optimization_result.at<gtsam::Point3>(Symbol('p', 4)));   
}
