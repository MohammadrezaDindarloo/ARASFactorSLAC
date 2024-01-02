#include "../include/libscampi.h"
#include <gtsam/slam/PriorFactor.h>
#include <gtsam/slam/BetweenFactor.h>

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
void forward_kinematic_factor_graph_optimizer(std::vector<Eigen::Matrix<double, 4, 1>> cable_length_collection,
                              std::vector<Eigen::Matrix<double, 2, 1>> cable_forces_collection , RobotParameters<double> params, 
                              std::vector<std::vector< Eigen::Matrix<double, 3, 1> >> b_in_w_collection, std::vector<std::vector< Eigen::Matrix<double, 3, 1> >> p_in_w_collection,
                              std::vector<Eigen::Matrix<double, 3, 1>> pulley_position_estimate,
                              std::vector<Eigen::Matrix<double, 3, 3>> rot_platform_collection,
                              gtsam::Rot3 DeltaRot, std::vector<Eigen::Matrix<double, 3, 1>> p_platform_collection, gtsam::Values *oprimization_result_LM)
{
    NonlinearFactorGraph graph;
    Values initial_estimate;

    // create noise for PTA mesurements
    // create a random engine with a seed
    std::default_random_engine generator(std::random_device{}());
    // create a uniform distribution between two number
    std::uniform_real_distribution<double> cable_1(-0.005, 0.005);
    std::uniform_real_distribution<double> cable_2(-0.005, 0.005);
    std::uniform_real_distribution<double> cable_3(-0.005, 0.005);
    std::uniform_real_distribution<double> cable_4(-0.005, 0.005);
    
    // Create a 6x6 matrix with 0.05 on the main diagonal and zeros elsewhere
    int rows = 6;
    int cols = 6;
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(rows, cols);
    matrix.diagonal() = 0.01 * Eigen::VectorXd::Ones(rows);

    auto Sensor_noiseModel_cost1 = gtsam::noiseModel::Isotropic::Sigma(4, sqrt(30));
    auto Sensor_noiseModel_cost2 = gtsam::noiseModel::Isotropic::Sigma(4, sqrt(500));
    auto Sensor_noiseModel_cost3 = gtsam::noiseModel::Isotropic::Sigma(4, sqrt(1));
    gtsam::noiseModel::Gaussian::shared_ptr noiseModel_pose3 = gtsam::noiseModel::Gaussian::Covariance(matrix);

    // Prior Factor for Robot Pose
    graph.emplace_shared<gtsam::PriorFactor<gtsam::Pose3>>(Symbol('X', 0), gtsam::Pose3(EigenMatrixToGtsamRot3(rot_platform_collection[0]), p_platform_collection[0]), noiseModel_pose3);
    
    for (size_t i = 0; i < p_platform_collection.size()-1; i++)
    {
        double random_noise_encoder_1 = cable_1(generator);
        double random_noise_encoder_2 = cable_2(generator);
        double random_noise_encoder_3 = cable_3(generator);
        double random_noise_encoder_4 = cable_4(generator);

        double encoder_data_1 = cable_length_collection[i][0] + random_noise_encoder_1 + 5.5;
        double encoder_data_2 = cable_length_collection[i][1] + random_noise_encoder_2 - 5.5;
        double encoder_data_3 = cable_length_collection[i][2] + random_noise_encoder_3 + 8.8;
        double encoder_data_4 = cable_length_collection[i][3] + random_noise_encoder_4 - 8.8;

        gtsam::Vector4 encodr_data = {encoder_data_1, encoder_data_2, encoder_data_3, encoder_data_4};

        graph.add(std::make_shared<FK_factor_graph_cost1>(Symbol('h', i), Symbol('v', i), Symbol('p', 1), Symbol('p', 2), Symbol('p', 3), Symbol('p', 4), Symbol('X', i), DeltaRot, Sensor_noiseModel_cost1));
        graph.add(std::make_shared<FK_factor_graph_cost2>(Symbol('h', i), Symbol('v', i), Symbol('o', 0), Symbol('o', 1), Symbol('o', 2), Symbol('o', 3), Symbol('p', 1), Symbol('p', 2), Symbol('p', 3), Symbol('p', 4), Symbol('X', i), encodr_data, DeltaRot, Sensor_noiseModel_cost2));
        graph.add(std::make_shared<FK_factor_graph_cost3>(Symbol('h', i), Symbol('v', i), Symbol('p', 1), Symbol('p', 2), Symbol('p', 3), Symbol('p', 4), Symbol('X', i), DeltaRot, Sensor_noiseModel_cost3));
        
        // Robot Pose Factors
        gtsam::Pose3 pose_k = gtsam::Pose3(EigenMatrixToGtsamRot3(rot_platform_collection[i]), p_platform_collection[i]);
        gtsam::Pose3 pose_k_plus1 = gtsam::Pose3(EigenMatrixToGtsamRot3(rot_platform_collection[i+1]), p_platform_collection[i+1]);
        gtsam::Pose3 relativePose = pose_k.between(pose_k_plus1);
        graph.emplace_shared<gtsam::BetweenFactor<gtsam::Pose3>>(Symbol('X', i), Symbol('X', i+1), relativePose, noiseModel_pose3);
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
    for (size_t i = 0; i < p_platform_collection.size()-1; i++)
    {
        initial_estimate.insert(Symbol('h', i), cable_forces_collection[i].norm()/sqrt(2));
        initial_estimate.insert(Symbol('v', i), cable_forces_collection[i].norm()/sqrt(2));  
        initial_estimate.insert(Symbol('X', i), gtsam::Pose3(EigenMatrixToGtsamRot3(rot_platform_collection[i]), p_platform_collection[i]));
    }
    initial_estimate.insert(Symbol('X', p_platform_collection.size()-1), gtsam::Pose3(EigenMatrixToGtsamRot3(rot_platform_collection[p_platform_collection.size()-1]), p_platform_collection[p_platform_collection.size()-1]));

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
void fkSolver(std::vector<Eigen::Matrix<double, 4, 1>> cable_length_collection,
              std::vector<Eigen::Matrix<double, 3, 1>> p_platform_collection,  
              std::vector<Eigen::Matrix<double, 3, 3>> rot_platform_collection, 
              std::vector<Eigen::Matrix<double, 3, 1>> pulley_position_estimate,
              std::vector<Eigen::Matrix<double, 2, 1>> cable_forces_collection, 
              RobotParameters<double> params, 
              FKDataOut<double> *result)
{   
    RobotState<double> state_init;
    GeometricVariables<double> geom_vars_calibration;
    std::vector<std::vector< Eigen::Matrix<double, 3, 1> >> b_in_w_collection;
    std::vector<std::vector< Eigen::Matrix<double, 3, 1> >> p_in_w_collection;
    for (size_t j = 0; j < p_platform_collection.size(); j++)
    {
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
    forward_kinematic_factor_graph_optimizer(cable_length_collection, cable_forces_collection, params, b_in_w_collection, p_in_w_collection,
                                             pulley_position_estimate, rot_platform_collection, 
                                             init_estimate_rotation, p_platform_collection, &optimization_result);

    // optimization_result.print();
    double disparity_in_f_h = 0.0;
    double disparity_in_f_v = 0.0;
    for (size_t i = 0; i < cable_forces_collection.size()-1; i++)
    {
        double temp1 = abs(optimization_result.at<double>(Symbol('h', i)) - cable_forces_collection[i][0]);
        double temp2 = abs(optimization_result.at<double>(Symbol('v', i)) - cable_forces_collection[i][1]);
        disparity_in_f_h = disparity_in_f_h + temp1;
        disparity_in_f_v = disparity_in_f_v + temp2;
    }
    std::cout << std::endl << "disparity_in_f_h: " << disparity_in_f_h / (cable_forces_collection.size()-1) << std::endl;
    std::cout << std::endl << "disparity_in_f_v: " << disparity_in_f_v / (cable_forces_collection.size()-1) << std::endl;
    //harvest the results
    result->pulley_position.push_back(optimization_result.at<gtsam::Point3>(Symbol('p', 1)));
    result->pulley_position.push_back(optimization_result.at<gtsam::Point3>(Symbol('p', 2)));
    result->pulley_position.push_back(optimization_result.at<gtsam::Point3>(Symbol('p', 3)));
    result->pulley_position.push_back(optimization_result.at<gtsam::Point3>(Symbol('p', 4)));   
}
