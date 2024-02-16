#include "../include/libscampi.h"
#include <gtsam/slam/PriorFactor.h>

// ****************************************** IK optimization *******************************************
// a function to create ivnerse kinematic factor graph 
void inverse_kinematic_factor_graph_optimizer(Eigen::Vector3d p_init, Eigen::Matrix3d rot_init, int largest_cable,
                              double init_estimate_h1, double init_estimate_v1, gtsam::Rot3 init_estimate_rot, gtsam::Values *oprimization_result_LM)
{
    NonlinearFactorGraph graph;
    Values initial_estimate;

    auto Sensor_noiseModel_cost1 = gtsam::noiseModel::Isotropic::Sigma(4, 0.001/3.0); // 0.0031/3
    auto Sensor_noiseModel_cost2 = gtsam::noiseModel::Isotropic::Sigma(4, 1.5/3.0); // 1.5/3
    auto Sensor_noiseModel_cost3 = gtsam::noiseModel::Isotropic::Sigma(4, 10.0); // 10.0

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
                                            std::vector<Eigen::Matrix<double, 4, 1>> IK_cable_length_collection,
                                            std::vector<Eigen::Matrix<double, 2, 1>> cable_forces_collection,
                                            std::vector<Eigen::Matrix<double, 3, 1>> p_platform_collection,
                                            std::vector<Eigen::Matrix<double, 3, 3>> rot_init_platform_collection,
                                            std::vector<Eigen::Matrix<double, 3, 3>> delta_rot_platform_collection,
                                            Eigen::Matrix<double, 4, 3> pulley_position_estimate,
                                            std::vector<std::vector< Eigen::Matrix<double, 3, 1> >> b_in_w_collection,
                                            std::vector<std::vector< Eigen::Matrix<double, 3, 1> >> p_in_w_collection,
                                            double *optimizer_error,
                                            std::vector<gtsam::Pose3> *Optimized_pose,
                                            std::vector<gtsam::Pose3> *GT_pose,
                                            gtsam::Values *oprimization_result_LM)
{
    NonlinearFactorGraph graph;
    Values initial_estimate;

    std::default_random_engine generator(std::random_device{}());

    // Encoder noise model
    double encoder_noise_gain = 0.0;
    std::normal_distribution<double> encoder_noise(0.0, encoder_noise_gain * (0.001)/3.0); // 0.1 degree in orientation error

    // UWB noise model
    double uwb_noise_gain = 0.0;
    std::normal_distribution<double> uwb_noise(0.0, uwb_noise_gain * (0.020)/3.0); // 0.1 degree in orientation error

    // Odomerty error
    gtsam::Point3 PositionNoise;
    double odometry_noise_gain = 0.0;
    std::normal_distribution<double> poistion_noise(0.0, odometry_noise_gain * (0.1/sqrt(3.0))/3.0); // 1 mm in position error
    std::normal_distribution<double> orientation_noise(0.0, odometry_noise_gain * (0.1/sqrt(3.0) * M_PI / 180.0)/3.0); // 0.1 degree in orientation error

    // Noise model for odometry and prior pose
    double translationnoise = 0.001/3.0; //in meter 
    double orientationnoise = (0.05 * M_PI / 180.0)/3.0; // in degree and convert to radian
    auto noiseModel_pose3 = noiseModel::Diagonal::Sigmas((gtsam::Vector(6)<<translationnoise,translationnoise,translationnoise,orientationnoise,orientationnoise,orientationnoise).finished());
    auto prior_noiseModel_pose3 = noiseModel::Diagonal::Sigmas((gtsam::Vector(6)<<0.05/3.0, 0.01/3.0, 0.003/3.0, (1.0 * M_PI / 180.0)/3.0, (5.0 * M_PI / 180.0)/3.0, (1.0 * M_PI / 180.0)/3.0).finished());
    auto prior_noiseModel_pose3_init = noiseModel::Diagonal::Sigmas((gtsam::Vector(6)<<1e-4, 1e-4, 1e-4, 1e-3, 1e-3, 1e-3).finished());

    // Cost noise models
    auto Sensor_noiseModel_cost1 = gtsam::noiseModel::Isotropic::Sigma(4, 0.001/3.0); // z    0.001    
    auto Sensor_noiseModel_cost2 = gtsam::noiseModel::Isotropic::Sigma(4, 0.001/3.0); // UWB       
    auto Sensor_noiseModel_cost3 = gtsam::noiseModel::Isotropic::Sigma(4, 0.005/3.0); // encoder   0.001
    
    auto pulley_location_noise_model = gtsam::noiseModel::Isotropic::Sigma(3, 1.0e-6); // z    0.001    

    std::vector<gtsam::Pose3> Optimized_pose_;
    std::vector<gtsam::Pose3> GT_pose_;
    int j = 0;
    int temp = 0;
    int Variaty = 0;
    std::vector<int> static_indexes;    
    for (size_t i = 0; i < p_platform_collection.size(); i++)
    {
        if (i>4 && std::abs(cable_length_collection[i][0]-cable_length_collection[i-1][0])<1.0e-3
                && std::abs(cable_length_collection[i][1]-cable_length_collection[i-1][1])<1.0e-3
                && std::abs(cable_length_collection[i][2]-cable_length_collection[i-1][2])<1.0e-3
                && std::abs(cable_length_collection[i][3]-cable_length_collection[i-1][3])<1.0e-3
                && std::abs(cable_length_collection[i][0]-cable_length_collection[i-2][0])<1.0e-3
                && std::abs(cable_length_collection[i][1]-cable_length_collection[i-2][1])<1.0e-3
                && std::abs(cable_length_collection[i][2]-cable_length_collection[i-2][2])<1.0e-3
                && std::abs(cable_length_collection[i][3]-cable_length_collection[i-2][3])<1.0e-3
                && std::abs(cable_length_collection[i][0]-cable_length_collection[i-3][0])<1.0e-3
                && std::abs(cable_length_collection[i][1]-cable_length_collection[i-3][1])<1.0e-3
                && std::abs(cable_length_collection[i][2]-cable_length_collection[i-3][2])<1.0e-3
                && std::abs(cable_length_collection[i][3]-cable_length_collection[i-3][3])<1.0e-3
                && std::abs((p_platform_collection[i]-p_platform_collection[i-1]).norm())<30.0e-3
                && std::abs((p_platform_collection[i]-p_platform_collection[i-2]).norm())<30.0e-3
                && std::abs((p_platform_collection[i]-p_platform_collection[i-3]).norm())<30.0e-3
                && p_platform_collection[i][0] > -0.10
                && p_platform_collection[i][0] <  0.71
                && p_platform_collection[i][1] > -3.14
                && p_platform_collection[i][1] <  0.54
                && p_platform_collection[i][2] <  4.20
                && p_platform_collection[i][2] >  0.00)
        {
            static_indexes.push_back(i);
        }
    }
    for (size_t i = 0; i < static_indexes.size(); i++)
    {
        std::cout << "static index: " << static_indexes[i] << std::endl;
    }

    graph.emplace_shared<gtsam::PriorFactor<gtsam::Pose3>>(Symbol('X', 0), gtsam::Pose3(EigenMatrixToGtsamRot3(rot_init_platform_collection[0]), p_platform_collection[0]), prior_noiseModel_pose3);
    for (int i = 0; i < p_platform_collection.size()-1; i++)
    {
        if (std::find(static_indexes.begin(), static_indexes.end(), i) != static_indexes.end() && 
            !(std::find(static_indexes.begin(), static_indexes.end(), i+1) != static_indexes.end()) 
            && std::abs(cable_length_collection[i][0]-IK_cable_length_collection[i][0]) < 0.005 &&
            std::abs(cable_length_collection[i][1]-IK_cable_length_collection[i][1]) < 0.005 &&
            std::abs(cable_length_collection[i][2]-IK_cable_length_collection[i][2]) < 0.005 &&
            std::abs(cable_length_collection[i][3]-IK_cable_length_collection[i][3]) < 0.005)
        {   
            std::cout << "static index: " << i << std::endl;
            std::cout << "diff_enc_inv: " << std::endl << cable_length_collection[i]-IK_cable_length_collection[i] << std::endl;
            if (i-temp != 1)
            {
                Variaty +=1;
            }
            temp = i;
            double uwb_data_1 = (b_in_w_collection[i][0] - p_in_w_collection[i][0]).norm() + uwb_noise(generator);
            double uwb_data_2 = (b_in_w_collection[i][1] - p_in_w_collection[i][1]).norm() + uwb_noise(generator);
            double uwb_data_3 = (b_in_w_collection[i][2] - p_in_w_collection[i][2]).norm() + uwb_noise(generator);
            double uwb_data_4 = (b_in_w_collection[i][3] - p_in_w_collection[i][3]).norm() + uwb_noise(generator);
            gtsam::Vector4 uwb_data = {uwb_data_1, uwb_data_2, uwb_data_3, uwb_data_4};

            double enc_data_1 = cable_length_collection[i][0];
            double enc_data_2 = cable_length_collection[i][1];
            double enc_data_3 = cable_length_collection[i][2];
            double enc_data_4 = cable_length_collection[i][3];
            gtsam::Vector4 enc_data = {enc_data_1, enc_data_2, enc_data_3, enc_data_4};

            graph.add(std::make_shared<FK_factor_graoh_cost1>(Symbol('h', j), Symbol('v', j), Symbol('r', j), Symbol('X', i), Symbol('p', 0), Symbol('p', 1), Symbol('p', 2), Symbol('p', 3), enc_data, Sensor_noiseModel_cost1));
            // graph.add(std::make_shared<FK_factor_graoh_cost2>(Symbol('h', j), Symbol('v', j), Symbol('r', j), Symbol('X', i), Symbol('p', 0), Symbol('p', 1), Symbol('p', 2), Symbol('p', 3), uwb_data, Sensor_noiseModel_cost2));
            graph.add(std::make_shared<FK_factor_graoh_cost3>(Symbol('h', j), Symbol('v', j), Symbol('r', j), Symbol('X', i), Symbol('p', 0), Symbol('p', 1), Symbol('p', 2), Symbol('p', 3), enc_data, Sensor_noiseModel_cost3));
       
            initial_estimate.insert(Symbol('h', j), cable_forces_collection[i][0]);
            initial_estimate.insert(Symbol('v', j), cable_forces_collection[i][1]);
            initial_estimate.insert(Symbol('r', j), EigenMatrixToGtsamRot3(delta_rot_platform_collection[i]));

            j +=1; 
        }

        // Robot Pose Factor
        // prior factor
        // if (i == 0)
        // {
        //     graph.emplace_shared<gtsam::PriorFactor<gtsam::Pose3>>(Symbol('X', i), gtsam::Pose3(EigenMatrixToGtsamRot3(rot_init_platform_collection[i]), p_platform_collection[i]), prior_noiseModel_pose3_init);
        // }
        // else
        // {
        // graph.emplace_shared<gtsam::PriorFactor<gtsam::Pose3>>(Symbol('X', i), gtsam::Pose3(EigenMatrixToGtsamRot3(rot_init_platform_collection[i]), p_platform_collection[i]), prior_noiseModel_pose3);
        // }
        // odometry factor
        gtsam::Pose3 pose_k = gtsam::Pose3(EigenMatrixToGtsamRot3(rot_init_platform_collection[i]), p_platform_collection[i]);
        GT_pose_.push_back(pose_k);
        gtsam::Pose3 pose_k_plus1 = gtsam::Pose3(EigenMatrixToGtsamRot3(rot_init_platform_collection[i+1]), p_platform_collection[i+1]);
        gtsam::Pose3 relativePose = pose_k.between(pose_k_plus1);
        // Generate random position noise
        PositionNoise.x() = poistion_noise(generator);
        PositionNoise.y() = poistion_noise(generator);
        PositionNoise.z() = poistion_noise(generator);
        gtsam::Point3 noisyPosition = relativePose.translation() + PositionNoise;
        // Generate random orientation noise
        gtsam::Rot3 orientationNoise = gtsam::Rot3::Expmap(gtsam::Vector3(orientation_noise(generator), orientation_noise(generator), orientation_noise(generator)));
        gtsam::Rot3 noisyOrientation = relativePose.rotation() * orientationNoise; 
        gtsam::Pose3 noisyRelativePose = gtsam::Pose3(noisyOrientation, noisyPosition);
        graph.emplace_shared<gtsam::BetweenFactor<gtsam::Pose3>>(Symbol('X', i), Symbol('X', i+1), noisyRelativePose, noiseModel_pose3);

        initial_estimate.insert(Symbol('X', i), gtsam::Pose3(EigenMatrixToGtsamRot3(rot_init_platform_collection[i]), p_platform_collection[i]));
    }

    initial_estimate.insert(Symbol('X', p_platform_collection.size()-1), gtsam::Pose3(EigenMatrixToGtsamRot3(rot_init_platform_collection[p_platform_collection.size()-1]), p_platform_collection[p_platform_collection.size()-1]));
    
    // graph.emplace_shared<gtsam::PriorFactor<gtsam::Point3>>(Symbol('p', 0), pulley_position_estimate.row(0), pulley_location_noise_model);
    // graph.emplace_shared<gtsam::PriorFactor<gtsam::Point3>>(Symbol('p', 1), pulley_position_estimate.row(1), pulley_location_noise_model);
    // graph.emplace_shared<gtsam::PriorFactor<gtsam::Point3>>(Symbol('p', 2), pulley_position_estimate.row(2), pulley_location_noise_model);
    // graph.emplace_shared<gtsam::PriorFactor<gtsam::Point3>>(Symbol('p', 3), pulley_position_estimate.row(3), pulley_location_noise_model);

    initial_estimate.insert(Symbol('p', 0), gtsam::Point3(pulley_position_estimate.row(0)));
    initial_estimate.insert(Symbol('p', 1), gtsam::Point3(pulley_position_estimate.row(1)));
    initial_estimate.insert(Symbol('p', 2), gtsam::Point3(pulley_position_estimate.row(2)));
    initial_estimate.insert(Symbol('p', 3), gtsam::Point3(pulley_position_estimate.row(3)));

    std::cout << "Number of static points: " << j << std::endl;
    std::cout << "Variaty of static points: " << Variaty << std::endl;

    gtsam::LevenbergMarquardtParams params; 
    int max_iterations = params.maxIterations;
    std::cout << std::endl;
    // params.setVerbosityLM("SUMMARY");
    LevenbergMarquardtOptimizer optimizer(graph, initial_estimate, params);
    Values result_LM = optimizer.optimize();
    for (size_t i = 0; i < p_platform_collection.size()-1; i++)
    {
        Optimized_pose_.push_back(result_LM.at<gtsam::Pose3>(Symbol('X', i)));
    }
    std::cout << std::endl << "forward kinematic optimization error: " << optimizer.error() << std::endl;
    *optimizer_error = optimizer.error();
    *oprimization_result_LM = result_LM;
    *Optimized_pose = Optimized_pose_;
    *GT_pose = GT_pose_;
}

// a function to hold parameters and invoke optimizer and back the optimized data
void fkSolver(  RobotParameters<double> params, 
                std::vector<Eigen::Matrix<double, 4, 1>> cable_length_collection, 
                std::vector<Eigen::Matrix<double, 4, 1>> IK_cable_length_collection,
                std::vector<Eigen::Matrix<double, 2, 1>> cable_forces_collection, 
                std::vector<Eigen::Matrix<double, 3, 1>> p_platform_collection,
                std::vector<Eigen::Matrix<double, 3, 3>> rot_init_platform_collection, 
                std::vector<Eigen::Matrix<double, 3, 3>> delta_rot_platform_collection,
                Eigen::Matrix<double, 4, 3> pulley_position_estimate,
                std::vector<gtsam::Pose3> *Optimized_pose,
                std::vector<gtsam::Pose3> *GT_pose,
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
    std::vector<gtsam::Pose3> Optimized_pose_;
    std::vector<gtsam::Pose3> GT_pose_;
    forward_kinematic_factor_graph_optimizer(cable_length_collection, IK_cable_length_collection, cable_forces_collection, p_platform_collection, 
                                             rot_init_platform_collection, delta_rot_platform_collection, pulley_position_estimate,
                                             b_in_w_collection, p_in_w_collection,
                                             &optimizer_error,
                                             &Optimized_pose_,
                                             &GT_pose_,
                                             &optimization_result);
    // optimization_result.print();
    //harvest the results
    Eigen::Matrix<double, 4, 3> pulley_optimized;
    *Optimized_pose = Optimized_pose_;
    *GT_pose = GT_pose_;

    pulley_optimized.row(0) = optimization_result.at<gtsam::Point3>(Symbol('p', 0));
    pulley_optimized.row(1) = optimization_result.at<gtsam::Point3>(Symbol('p', 1));
    pulley_optimized.row(2) = optimization_result.at<gtsam::Point3>(Symbol('p', 2));
    pulley_optimized.row(3) = optimization_result.at<gtsam::Point3>(Symbol('p', 3));

    result->pulley_result = pulley_optimized;
    result->optimizer_error = optimizer_error; 
}
