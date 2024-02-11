#include "../include/libscampi.h"
#include <gtsam/slam/PriorFactor.h>

// ****************************************** IK optimization *******************************************
// a function to create ivnerse kinematic factor graph 
void inverse_kinematic_factor_graph_optimizer(Eigen::Vector3d p_init, Eigen::Matrix3d rot_init, int largest_cable,
                              double init_estimate_h1, double init_estimate_v1, gtsam::Rot3 init_estimate_rot, gtsam::Values *oprimization_result_LM)
{
    NonlinearFactorGraph graph;
    Values initial_estimate;

    auto Sensor_noiseModel_cost1 = gtsam::noiseModel::Isotropic::Sigma(4, 0.002);
    auto Sensor_noiseModel_cost2 = gtsam::noiseModel::Isotropic::Sigma(4, 10.0);
    auto Sensor_noiseModel_cost3 = gtsam::noiseModel::Isotropic::Sigma(4, 100.0);

    graph.add(std::make_shared<IK_factor_graoh_cost1>(Symbol('h', 1), Symbol('v', 1), Symbol('r', 1), p_init, rot_init, largest_cable, Sensor_noiseModel_cost1));
    graph.add(std::make_shared<IK_factor_graoh_cost2>(Symbol('h', 1), Symbol('v', 1), Symbol('r', 1), p_init, rot_init, largest_cable, Sensor_noiseModel_cost2));
    graph.add(std::make_shared<IK_factor_graoh_cost3>(Symbol('h', 1), Symbol('v', 1), Symbol('r', 1), p_init, rot_init, largest_cable, Sensor_noiseModel_cost3));

    initial_estimate.insert(Symbol('h', 1), init_estimate_h1);
    initial_estimate.insert(Symbol('v', 1), init_estimate_v1);
    initial_estimate.insert(Symbol('r', 1), init_estimate_rot);

    gtsam::LevenbergMarquardtParams params; 
    int max_iterations = params.maxIterations;
    LevenbergMarquardtOptimizer optimizer(graph, initial_estimate, params);
    Values result_LM = optimizer.optimize();

    std::cout << std::endl << "inverse kinematic optimization error: " << optimizer.error() << std::endl;
    if (optimizer.error() > 1e0)
    {
        std::cout << std::endl << "!!!!!!!!!!!!!! inverse kinematic error is too high !!!!!!!!!!!!!!" << std::endl;
        std::cout << std::endl << "!!!!!!!!!!!!!! p_init !!!!!!!!!!!!!!" << p_init << std::endl;
    }
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
    inverse_kinematic_factor_graph_optimizer(p_platform, rot_init, largest_cable,
                                            init_estimate_h1, init_estimate_v1, init_estimate_rot, &optimization_result);

    // optimization_result.print();
    //harvest the results
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
                                            std::vector<Eigen::Matrix<double, 2, 1>> cable_forces_collection,
                                            std::vector<Eigen::Matrix<double, 3, 1>> p_platform_collection,
                                            std::vector<Eigen::Matrix<double, 3, 3>> rot_init_platform_collection,
                                            std::vector<Eigen::Matrix<double, 3, 3>> delta_rot_platform_collection,
                                            Eigen::Matrix<double, 4, 3> pulley_position_estimate,
                                            std::vector<std::vector< Eigen::Matrix<double, 3, 1> >> b_in_w_collection,
                                            std::vector<std::vector< Eigen::Matrix<double, 3, 1> >> p_in_w_collection,
                                            double *optimizer_error,
                                            gtsam::Values *oprimization_result_LM)
{
    NonlinearFactorGraph graph;
    Values initial_estimate;
    double pulley_noise = 1e-5;
    // auto Sensor_noiseModel_cost1 = gtsam::noiseModel::Isotropic::Sigma(4, 0.001); // z          0.001
    // auto Sensor_noiseModel_cost2 = gtsam::noiseModel::Isotropic::Sigma(4, 0.002); // UWB        0.001
    // auto Sensor_noiseModel_cost3 = gtsam::noiseModel::Isotropic::Sigma(4, 0.0004); // encoder   0.0002

    auto Sensor_noiseModel_cost1 = noiseModel::Diagonal::Sigmas((gtsam::Vector(4)<<0.0003, 0.0003, 0.0003, 0.0003).finished());      // z 
    auto Sensor_noiseModel_cost2 = noiseModel::Diagonal::Sigmas((gtsam::Vector(4)<<0.001/3.0, 0.001/3.0, 0.001/3.0, 0.001/3.0).finished());      // UWB
    auto Sensor_noiseModel_cost3 = noiseModel::Diagonal::Sigmas((gtsam::Vector(4)<<0.001/3.0, 0.001/3.0, 0.001/3.0, 0.001/3.0).finished());  // encoder

    auto prior_noiseModel_pulley = noiseModel::Diagonal::Sigmas((gtsam::Vector(3)<<pulley_noise, pulley_noise, pulley_noise).finished());

    for (size_t i = 0; i < p_platform_collection.size(); i++)
    {
        double uwb_data_1 = (b_in_w_collection[i][0] - p_in_w_collection[i][0]).norm();
        double uwb_data_2 = (b_in_w_collection[i][1] - p_in_w_collection[i][1]).norm();
        double uwb_data_3 = (b_in_w_collection[i][2] - p_in_w_collection[i][2]).norm();
        double uwb_data_4 = (b_in_w_collection[i][3] - p_in_w_collection[i][3]).norm();
        gtsam::Vector4 uwb_data = {uwb_data_1, uwb_data_2, uwb_data_3, uwb_data_4};

        graph.add(std::make_shared<FK_factor_graoh_cost1>(Symbol('h', i), Symbol('v', i), Symbol('r', i), Symbol('t', i), Symbol('p', 0), Symbol('p', 1), Symbol('p', 2), Symbol('p', 3), cable_length_collection[i], rot_init_platform_collection[i], Sensor_noiseModel_cost1));
        // graph.add(std::make_shared<FK_factor_graoh_cost2>(Symbol('h', i), Symbol('v', i), Symbol('r', i), Symbol('t', i), Symbol('p', 0), Symbol('p', 1), Symbol('p', 2), Symbol('p', 3), uwb_data, rot_init_platform_collection[i], Sensor_noiseModel_cost2));
        graph.add(std::make_shared<FK_factor_graoh_cost3>(Symbol('h', i), Symbol('v', i), Symbol('r', i), Symbol('t', i), Symbol('p', 0), Symbol('p', 1), Symbol('p', 2), Symbol('p', 3), cable_length_collection[i], rot_init_platform_collection[i], Sensor_noiseModel_cost3));

        initial_estimate.insert(Symbol('h', i), cable_forces_collection[i][0]);
        initial_estimate.insert(Symbol('v', i), cable_forces_collection[i][1]);
        initial_estimate.insert(Symbol('r', i), EigenMatrixToGtsamRot3(delta_rot_platform_collection[i]));
        initial_estimate.insert(Symbol('t', i), p_platform_collection[i]);
    }

    // graph.add(gtsam::PriorFactor<gtsam::Point3>(Symbol('p', 0), pulley_position_estimate.row(0), prior_noiseModel_pulley));
    // graph.add(gtsam::PriorFactor<gtsam::Point3>(Symbol('p', 1), pulley_position_estimate.row(1), prior_noiseModel_pulley));
    // graph.add(gtsam::PriorFactor<gtsam::Point3>(Symbol('p', 2), pulley_position_estimate.row(2), prior_noiseModel_pulley));
    // graph.add(gtsam::PriorFactor<gtsam::Point3>(Symbol('p', 3), pulley_position_estimate.row(3), prior_noiseModel_pulley));

    initial_estimate.insert(Symbol('p', 0), gtsam::Point3(pulley_position_estimate.row(0)));
    initial_estimate.insert(Symbol('p', 1), gtsam::Point3(pulley_position_estimate.row(1)));
    initial_estimate.insert(Symbol('p', 2), gtsam::Point3(pulley_position_estimate.row(2)));
    initial_estimate.insert(Symbol('p', 3), gtsam::Point3(pulley_position_estimate.row(3)));

    gtsam::LevenbergMarquardtParams params; 
    params.relativeErrorTol = 1e-23;
    int max_iterations = params.maxIterations;
    LevenbergMarquardtOptimizer optimizer(graph, initial_estimate, params);
    Values result_LM = optimizer.optimize();

    // // Set up variables for iteration and convergence check
    // double convergence_threshold = 1e-9; // Set your desired convergence threshold
    // int iteration = 0;
    // while (iteration < max_iterations) 
    // {
    //     // Run one iteration of the optimizer
    //     optimizer.iterate();
    //     // Get the current error after the iteration
    //     double current_error = optimizer.error();
    //     // Print the error for this iteration
    //     std::cout << "Iteration " << iteration + 1 << ": Error = " << current_error << std::endl;
    //     // Check for convergence
    //     if (current_error < convergence_threshold) {
    //         std::cout << "Converged at iteration " << iteration + 1 << std::endl;
    //         break;
    //     }
    //     iteration++;
    // }

    std::cout << std::endl << "forward kinematic optimization error: " << optimizer.error() << std::endl;
    *optimizer_error = optimizer.error();
    *oprimization_result_LM = result_LM;
}

// a function to hold parameters and invoke optimizer and back the optimized data
void fkSolver(  RobotParameters<double> params, 
                std::vector<Eigen::Matrix<double, 4, 1>> cable_length_collection, 
                std::vector<Eigen::Matrix<double, 2, 1>> cable_forces_collection, 
                std::vector<Eigen::Matrix<double, 3, 1>> p_platform_collection,
                std::vector<Eigen::Matrix<double, 3, 3>> rot_init_platform_collection, 
                std::vector<Eigen::Matrix<double, 3, 3>> delta_rot_platform_collection,
                Eigen::Matrix<double, 4, 3> pulley_position_estimate,
                FKDataOut<double> *result)
{  
    RobotState<double> state_init;
    GeometricVariables<double> geom_vars_calibration;
    std::vector<std::vector< Eigen::Matrix<double, 3, 1> >> b_in_w_collection;
    std::vector<std::vector< Eigen::Matrix<double, 3, 1> >> p_in_w_collection;
    for (size_t j = 0; j < p_platform_collection.size(); j++)
    {
        state_init.p_platform = p_platform_collection[j];
        state_init.rot_platform = rot_init_platform_collection[j] * delta_rot_platform_collection[j];
        getGeometricVariables<double>(state_init, params, &geom_vars_calibration);
        b_in_w_collection.push_back(geom_vars_calibration.b_in_w);
        p_in_w_collection.push_back(geom_vars_calibration.p_in_w);
    }

    // run optimization!
    double optimizer_error;
    gtsam::Values optimization_result; 
    forward_kinematic_factor_graph_optimizer(cable_length_collection, cable_forces_collection, p_platform_collection, 
                                             rot_init_platform_collection, delta_rot_platform_collection, pulley_position_estimate,
                                             b_in_w_collection, p_in_w_collection,
                                             &optimizer_error,
                                             &optimization_result);
    // optimization_result.print();
    //harvest the results
    Eigen::Matrix<double, 4, 3> pulley_optimized;
    
    pulley_optimized.row(0) = optimization_result.at<gtsam::Point3>(Symbol('p', 0));
    pulley_optimized.row(1) = optimization_result.at<gtsam::Point3>(Symbol('p', 1));
    pulley_optimized.row(2) = optimization_result.at<gtsam::Point3>(Symbol('p', 2));
    pulley_optimized.row(3) = optimization_result.at<gtsam::Point3>(Symbol('p', 3));

    Eigen::Matrix3d delta_rot_optimized = gtsamRot3ToEigenMatrix(optimization_result.at<gtsam::Rot3>(Symbol('r', 0)));
    Matrix3d rot_result = rot_init_platform_collection[0] * delta_rot_optimized;
    result->rot_platform = rot_result;
 
    GeometricVariables<double> geom_vars;
    CatenaryVariables<double> cat_vars;
    RobotState<double> state;
    state.rot_platform = rot_result;
    Eigen::VectorXd p_platform = optimization_result.at<gtsam::Point3>(Symbol('t', 0));
    result->p_platform = p_platform;
    state.p_platform = result->p_platform;

    double fh = optimization_result.at<double>(Symbol('h', 0));
    double fv = optimization_result.at<double>(Symbol('v', 0));

    getGeometricVariables<double>(state,params,&geom_vars);
    getCableForces<double>(fh, fv, &state, params,geom_vars);
    getCatenaryVariables<double>(state,params, geom_vars,&cat_vars);

    result->c1 = cat_vars.c1;
    result->c2 = cat_vars.c2;
    result->lc_cat = cat_vars.lc_cat;
    result->cable_forces = state.cable_forces;
    result->b_in_w = geom_vars.b_in_w;
    result->pulley_result = pulley_optimized;
    result->optimizer_error = optimizer_error; 
}
