#include "../include/libscampi.h"
#include <gtsam/slam/PriorFactor.h>

// ****************************************** IK optimization *******************************************
// a function to create ivnerse kinematic factor graph 
void inverse_kinematic_factor_graph_optimizer(Eigen::Vector3d p_init, Eigen::Matrix3d rot_init, int largest_cable,
                              double init_estimate_h1, double init_estimate_v1, gtsam::Rot3 init_estimate_rot, Eigen::Matrix<double, 4, 3> pulley_position_estimate, 
                              int inner_interval, Eigen::Matrix<double, 4, 1> cable_length, VectorXi reorder_idx, gtsam::Values *oprimization_result_LM)
{
    Eigen::Matrix<double, 4, 1> cable_length_reorder;
    for(int i = 0; i < cable_length.size(); i++)
    {
        cable_length_reorder[i] = cable_length[reorder_idx[i]];
    }

    NonlinearFactorGraph graph;
    Values initial_estimate;
    double interval_rate = 1.0; //std::pow((1.0/double(interval_rate + 1)),2)
    auto Sensor_noiseModel_cost1 = gtsam::noiseModel::Isotropic::Sigma(4, 0.205/3.0); // 0.02/3.0
    auto Sensor_noiseModel_cost2 = gtsam::noiseModel::Isotropic::Sigma(4, 0.5/3.0); // 0.02/3.0
    auto Sensor_noiseModel_cost3 = gtsam::noiseModel::Isotropic::Sigma(4, 1.0/3.0); // 1.0/3.0
    auto prior_noiseModel_delta_rot = noiseModel::Diagonal::Sigmas((gtsam::Vector(3)<< 30.0 * M_PI/180.0,  30.0 * M_PI/180.0,  30.0 * M_PI/180.0).finished()); // 2.0e-1
    auto prior_noiseModel_pulley = noiseModel::Diagonal::Sigmas((gtsam::Vector(3)<< 5.0/sqrt(3.0)/3.0,  5.0/sqrt(3.0)/3.0,  5.0/sqrt(3.0)/3.0).finished());

    graph.add(std::make_shared<IK_factor_graoh_cost1>(Symbol('h', 1), Symbol('v', 1), Symbol('r', 1), Symbol('p', 0), Symbol('p', 1), Symbol('p', 2), Symbol('p', 3), p_init, rot_init, largest_cable, cable_length, Sensor_noiseModel_cost1));
    graph.add(std::make_shared<IK_factor_graoh_cost2>(Symbol('h', 1), Symbol('v', 1), Symbol('r', 1), Symbol('p', 0), Symbol('p', 1), Symbol('p', 2), Symbol('p', 3), p_init, rot_init, largest_cable, cable_length, Sensor_noiseModel_cost2));
    graph.add(std::make_shared<IK_factor_graoh_cost3>(Symbol('h', 1), Symbol('v', 1), Symbol('r', 1), Symbol('p', 0), Symbol('p', 1), Symbol('p', 2), Symbol('p', 3), p_init, rot_init, largest_cable, cable_length, Sensor_noiseModel_cost3));

    // graph.emplace_shared<gtsam::PriorFactor<gtsam::Rot3>>(Symbol('r', 1), init_estimate_rot, prior_noiseModel_delta_rot);

    graph.emplace_shared<gtsam::PriorFactor<gtsam::Point3>>(Symbol('p', 0), gtsam::Point3(pulley_position_estimate.row(0)), prior_noiseModel_pulley);
    graph.emplace_shared<gtsam::PriorFactor<gtsam::Point3>>(Symbol('p', 1), gtsam::Point3(pulley_position_estimate.row(1)), prior_noiseModel_pulley);
    graph.emplace_shared<gtsam::PriorFactor<gtsam::Point3>>(Symbol('p', 2), gtsam::Point3(pulley_position_estimate.row(2)), prior_noiseModel_pulley);
    graph.emplace_shared<gtsam::PriorFactor<gtsam::Point3>>(Symbol('p', 3), gtsam::Point3(pulley_position_estimate.row(3)), prior_noiseModel_pulley);

    initial_estimate.insert(Symbol('h', 1), init_estimate_h1);
    initial_estimate.insert(Symbol('v', 1), init_estimate_v1);
    initial_estimate.insert(Symbol('r', 1), init_estimate_rot);
    
    initial_estimate.insert(Symbol('p', 0), gtsam::Point3(pulley_position_estimate.row(0)));
    initial_estimate.insert(Symbol('p', 1), gtsam::Point3(pulley_position_estimate.row(1)));
    initial_estimate.insert(Symbol('p', 2), gtsam::Point3(pulley_position_estimate.row(2)));
    initial_estimate.insert(Symbol('p', 3), gtsam::Point3(pulley_position_estimate.row(3)));

    gtsam::LevenbergMarquardtParams params; 
    int max_iterations = params.maxIterations;
    LevenbergMarquardtOptimizer optimizer(graph, initial_estimate, params);
    Values result_LM = optimizer.optimize();

    std::cout << "inverse kinematic optimization error: " << optimizer.error() << std::endl;
    std::cout << "inverse kinematic optimization error: " << optimizer.error() << std::endl;
    if (optimizer.error() > 1e4)
    {
        std::cout << "!!!!!!!!!!!!!! inverse kinematic error is too high !!!!!!!!!!!!!!" << std::endl;
    }
    *oprimization_result_LM = result_LM;
    
    // std::cout << "p0_inv: " << std::endl << result_LM.at<gtsam::Point3>(Symbol('p', 0)) << std::endl;
    // std::cout << "p1_inv: " << std::endl << result_LM.at<gtsam::Point3>(Symbol('p', 1)) << std::endl;
    // std::cout << "p2_inv: " << std::endl << result_LM.at<gtsam::Point3>(Symbol('p', 2)) << std::endl;
    // std::cout << "p3_inv: " << std::endl << result_LM.at<gtsam::Point3>(Symbol('p', 3)) << std::endl;
}

// a function to hold parameters and invoke optimizer and back the optimized data
void ikSolver(RobotParameters<double> params, 
              Eigen::Matrix3d rot_init, 
              Eigen::Vector3d p_platform,
              Eigen::Matrix<double, 4, 1> cable_length,
              Eigen::Matrix<double, 4, 3> pulley_position_estimate,
              int inner_interval,
              IKDataOut<double> *result)
{
    //reorder the cable forces and choose the cable with largest length as the the first cable (for numerical stability)
    VectorXi reorder_idx(params.pulleys.size());
    RobotParameters<double> params_reord;
    RobotState<double> state;
    state.rot_platform = rot_init;
    state.p_platform = p_platform;
    changeOrderForSolver<double>(state, params, cable_length, &params_reord, &reorder_idx);
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
    // gtsam::Rot3 init_estimate_rot = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    gtsam::Rot3 init_estimate_rot = gtsam::Rot3();
    // run optimization!
    gtsam::Values optimization_result;
    inverse_kinematic_factor_graph_optimizer(p_platform, rot_init, largest_cable,
                                            init_estimate_h1, init_estimate_v1, init_estimate_rot, pulley_position_estimate, inner_interval, cable_length, reorder_idx, &optimization_result);

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
void forward_kinematic_factor_graph_optimizer(std::vector<double> cable_offset,
                                            std::vector<Eigen::Matrix<double, 4, 1>> cable_length_collection,
                                            std::vector<Eigen::Matrix<double, 2, 1>> cable_forces_collection,
                                            std::vector<Eigen::Matrix<double, 3, 1>> p_platform_collection,
                                            std::vector<Eigen::Matrix<double, 3, 3>> rot_init_platform_collection,
                                            std::vector<Eigen::Matrix<double, 3, 3>> delta_rot_platform_collection,
                                            Eigen::Matrix<double, 4, 3> pulley_position_estimate,
                                            std::vector<std::vector< Eigen::Matrix<double, 3, 1> >> b_in_w_collection,
                                            std::vector<std::vector< Eigen::Matrix<double, 3, 1> >> p_in_w_collection,
                                            std::vector<double> first_cable_force_magnitude,
                                            double *optimizer_error,
                                            std::vector<gtsam::Pose3> *Optimized_pose,
                                            std::vector<gtsam::Pose3> *GT_pose,
                                            gtsam::Values *oprimization_result_LM)
{
    NonlinearFactorGraph graph;
    Values initial_estimate;
    auto sensor_mode_prior_factory = true;
    std::default_random_engine generator(std::random_device{}());

    // noise hander
    double all_noise_activation = 0.0;

    // Encoder noise model
    double encoder_noise_gain = 1.0;
    encoder_noise_gain = encoder_noise_gain * all_noise_activation;
    std::normal_distribution<double> encoder_noise(0.0, encoder_noise_gain * (0.01)/3.0); // 0.1 degree in orientation error

    // UWB noise model
    double uwb_noise_gain = 1.0;
    uwb_noise_gain = uwb_noise_gain * all_noise_activation;
    std::normal_distribution<double> uwb_noise(0.0, uwb_noise_gain * (0.020)/3.0); // 0.1 degree in orientation error

    gtsam::Point3 PositionNoise;
    gtsam::Rot3 orientationNoise;
    auto prior_noiseModel_init_pose = noiseModel::Diagonal::Sigmas((gtsam::Vector(6)<<1e-6, 1e-6, 1e-6, 1e-4, 1e-4, 1e-4).finished());

    // Noise model for odometry pose
    double odometry_noise_gain_odometry = 1.0;
    odometry_noise_gain_odometry = odometry_noise_gain_odometry * all_noise_activation;
    std::normal_distribution<double> poistion_noise_odometry(0.0, odometry_noise_gain_odometry * (0.001/sqrt(3.0))/3.0); // 1 mm in position error
    std::normal_distribution<double> orientation_noise_odometry(0.0, odometry_noise_gain_odometry * (0.1/sqrt(3.0) * M_PI / 180.0)/3.0); // 0.1 degree in orientation error
    double translationnoise_odometry = 0.001/sqrt(3.0)/3.0; //in meter 
    double orientationnoise_odometry = (0.1/sqrt(3.0) * M_PI / 180.0)/3.0; // in degree and convert to radian
    auto noiseModel_pose3 = noiseModel::Diagonal::Sigmas((gtsam::Vector(6)<<translationnoise_odometry,translationnoise_odometry,translationnoise_odometry,orientationnoise_odometry,orientationnoise_odometry,orientationnoise_odometry).finished());

    // Noise model for prior pose
    double odometry_noise_gain_prior = 1.0;
    odometry_noise_gain_prior = odometry_noise_gain_prior * all_noise_activation;
    std::normal_distribution<double> poistion_noise_prior(0.0, odometry_noise_gain_prior * (0.01/sqrt(3.0))/3.0); // 1 mm in position error
    std::normal_distribution<double> orientation_noise_prior(0.0, odometry_noise_gain_prior * (1.0/sqrt(3.0) * M_PI / 180.0)/3.0); // 0.1 degree in orientation error
    double translationnoise_prior = (0.001/sqrt(3.0))/3.0; //in meter 
    double orientationnoise_prior = (0.1/sqrt(3.0) * M_PI / 180.0)/3.0; // in degree and convert to radian
    auto prior_noiseModel_pose3 = noiseModel::Diagonal::Sigmas((gtsam::Vector(6)<<translationnoise_prior, translationnoise_prior, translationnoise_prior, orientationnoise_prior, orientationnoise_prior, orientationnoise_prior).finished());
    
    auto prior_noiseModel_offset = noiseModel::Diagonal::Sigmas((gtsam::Vector(4)<<1e-4, 1e-4, 1e-4, 1e-4).finished());
    auto prior_noiseModel_delta_rot = noiseModel::Diagonal::Sigmas((gtsam::Vector(3)<< (30.0e0 * M_PI/180.0) / 3.0, (30.0e0 * M_PI/180.0) / 3.0, (30.0e0 * M_PI/180.0) / 3.0).finished()); // 1.0e-2
    auto prior_noiseModel_pulley = noiseModel::Diagonal::Sigmas((gtsam::Vector(3)<< 5.0/sqrt(3.0)/3.0,  5.0/sqrt(3.0)/3.0,  5.0/sqrt(3.0)/3.0).finished());

    // Cost noise models
    auto Sensor_noiseModel_cost1 = noiseModel::Diagonal::Sigmas((gtsam::Vector(4)<< 0.05/3.0, 0.05/3.0, 0.05/3.0, 0.05/3.0).finished());
    // auto Sensor_noiseModel_cost1 = gtsam::noiseModel::Isotropic::Sigma(4, 0.15/3.0 ); // z     1.5
    auto Sensor_noiseModel_cost2 = gtsam::noiseModel::Isotropic::Sigma(4, 1.0/3.0); // l -||b-a||   
    auto Sensor_noiseModel_cost3 = gtsam::noiseModel::Isotropic::Sigma(4, 0.05/3.0); // encoder   0.35

    std::vector<gtsam::Pose3> Optimized_pose_;
    std::vector<gtsam::Pose3> GT_pose_;
    graph.emplace_shared<gtsam::PriorFactor<gtsam::Pose3>>(Symbol('X', 0), gtsam::Pose3(EigenMatrixToGtsamRot3(rot_init_platform_collection[0]), p_platform_collection[0]), prior_noiseModel_init_pose);
    for (size_t i = 0; i < p_platform_collection.size(); i++)
    {
        double enc_data_1 = cable_length_collection[i][0] + encoder_noise(generator);
        double enc_data_2 = cable_length_collection[i][1] + encoder_noise(generator);
        double enc_data_3 = cable_length_collection[i][2] + encoder_noise(generator);
        double enc_data_4 = cable_length_collection[i][3] + encoder_noise(generator);
        gtsam::Vector4 enc_data = {enc_data_1, enc_data_2, enc_data_3, enc_data_4};

        graph.add(std::make_shared<FK_factor_graoh_cost1>(Symbol('h', i), Symbol('v', i), Symbol('r', i), Symbol('X', i), Symbol('p', 0), Symbol('p', 1), Symbol('p', 2), Symbol('p', 3), enc_data, Sensor_noiseModel_cost1));
        // graph.add(std::make_shared<FK_factor_graoh_cost2>(Symbol('h', i), Symbol('v', i), Symbol('r', i), Symbol('X', i), Symbol('p', 0), Symbol('p', 1), Symbol('p', 2), Symbol('p', 3), enc_data, Sensor_noiseModel_cost2));
        graph.add(std::make_shared<FK_factor_graoh_cost3>(Symbol('h', i), Symbol('v', i), Symbol('r', i), Symbol('X', i), Symbol('p', 0), Symbol('p', 1), Symbol('p', 2), Symbol('p', 3), Symbol('o', 0), enc_data, Sensor_noiseModel_cost3));
        
        if (sensor_mode_prior_factory)
        {
            // ******** prior factor ********
            gtsam::Pose3 pose_k;
            gtsam::Pose3 pose_k_GT;
            if (i == 0)
            {
                pose_k = gtsam::Pose3(EigenMatrixToGtsamRot3(rot_init_platform_collection[i]), p_platform_collection[i]);
                pose_k_GT = pose_k;
            }
            else
            {
                // Generate random position noise
                PositionNoise.x() = poistion_noise_prior(generator);
                PositionNoise.y() = poistion_noise_prior(generator);
                PositionNoise.z() = poistion_noise_prior(generator);
                gtsam::Point3 noisyPosition = p_platform_collection[i] + PositionNoise;
                // Generate random orientation noise
                orientationNoise = gtsam::Rot3::Expmap(gtsam::Vector3(orientation_noise_prior(generator), orientation_noise_prior(generator), orientation_noise_prior(generator)));
                gtsam::Rot3 noisyOrientation = EigenMatrixToGtsamRot3(rot_init_platform_collection[i]) * orientationNoise;
                pose_k = gtsam::Pose3(noisyOrientation, noisyPosition);
                pose_k_GT = gtsam::Pose3(EigenMatrixToGtsamRot3(rot_init_platform_collection[i]), p_platform_collection[i]);
                graph.emplace_shared<gtsam::PriorFactor<gtsam::Pose3>>(Symbol('X', i), pose_k, prior_noiseModel_pose3);
            }
            GT_pose_.push_back(pose_k);
        }

        else
        {
            // ******** odometry factor ********
            gtsam::Pose3 pose_k = gtsam::Pose3(EigenMatrixToGtsamRot3(rot_init_platform_collection[i]), p_platform_collection[i]);
            GT_pose_.push_back(pose_k);
            gtsam::Pose3 pose_k_plus1 = gtsam::Pose3(EigenMatrixToGtsamRot3(rot_init_platform_collection[i+1]), p_platform_collection[i+1]);
            gtsam::Pose3 relativePose = pose_k.between(pose_k_plus1);
            // Generate random position noise
            PositionNoise.x() = poistion_noise_odometry(generator);
            PositionNoise.y() = poistion_noise_odometry(generator);
            PositionNoise.z() = poistion_noise_odometry(generator);
            gtsam::Point3 noisyPosition = relativePose.translation() + PositionNoise;
            // Generate random orientation noise
            orientationNoise = gtsam::Rot3::Expmap(gtsam::Vector3(orientation_noise_odometry(generator), orientation_noise_odometry(generator), orientation_noise_odometry(generator)));
            gtsam::Rot3 noisyOrientation = relativePose.rotation() * orientationNoise; 
            gtsam::Pose3 noisyRelativePose = gtsam::Pose3(noisyOrientation, noisyPosition);
            graph.emplace_shared<gtsam::BetweenFactor<gtsam::Pose3>>(Symbol('X', i), Symbol('X', i+1), noisyRelativePose, noiseModel_pose3);
            if (i==p_platform_collection.size()-1)
            {
                initial_estimate.insert(Symbol('X', i+1), gtsam::Pose3(EigenMatrixToGtsamRot3(rot_init_platform_collection[i+1]), p_platform_collection[i+1]));
            }
        }

        // graph.emplace_shared<gtsam::PriorFactor<gtsam::Rot3>>(Symbol('r', i), EigenMatrixToGtsamRot3(delta_rot_platform_collection[i]), prior_noiseModel_delta_rot);

        initial_estimate.insert(Symbol('h', i), cable_forces_collection[i][0]);
        initial_estimate.insert(Symbol('v', i), cable_forces_collection[i][1]);
        initial_estimate.insert(Symbol('r', i), EigenMatrixToGtsamRot3(delta_rot_platform_collection[i])); 
        initial_estimate.insert(Symbol('X', i), gtsam::Pose3(EigenMatrixToGtsamRot3(rot_init_platform_collection[i]), p_platform_collection[i]));
    }

    graph.emplace_shared<gtsam::PriorFactor<gtsam::Vector4>>(Symbol('o', 0), gtsam::Vector4(0.0, 0.0, 0.0, 0.0), prior_noiseModel_offset);
       
    // graph.emplace_shared<gtsam::PriorFactor<gtsam::Point3>>(Symbol('p', 0), gtsam::Point3(pulley_position_estimate.row(0)), prior_noiseModel_pulley);
    // graph.emplace_shared<gtsam::PriorFactor<gtsam::Point3>>(Symbol('p', 1), gtsam::Point3(pulley_position_estimate.row(1)), prior_noiseModel_pulley);
    // graph.emplace_shared<gtsam::PriorFactor<gtsam::Point3>>(Symbol('p', 2), gtsam::Point3(pulley_position_estimate.row(2)), prior_noiseModel_pulley);
    // graph.emplace_shared<gtsam::PriorFactor<gtsam::Point3>>(Symbol('p', 3), gtsam::Point3(pulley_position_estimate.row(3)), prior_noiseModel_pulley);

    // Initial estimate
    initial_estimate.insert(Symbol('p', 0), gtsam::Point3(pulley_position_estimate.row(0)));
    initial_estimate.insert(Symbol('p', 1), gtsam::Point3(pulley_position_estimate.row(1)));
    initial_estimate.insert(Symbol('p', 2), gtsam::Point3(pulley_position_estimate.row(2)));
    initial_estimate.insert(Symbol('p', 3), gtsam::Point3(pulley_position_estimate.row(3)));
    initial_estimate.insert(Symbol('o', 0), gtsam::Vector4(0.0, 0.0, 0.0, 0.0));

    gtsam::LevenbergMarquardtParams params; 
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
                std::vector<double> cable_offset,
                std::vector<Eigen::Matrix<double, 4, 1>> cable_length_collection, 
                std::vector<Eigen::Matrix<double, 2, 1>> cable_forces_collection, 
                std::vector<Eigen::Matrix<double, 3, 1>> p_platform_collection,
                std::vector<Eigen::Matrix<double, 3, 3>> rot_init_platform_collection, 
                std::vector<Eigen::Matrix<double, 3, 3>> delta_rot_platform_collection,
                Eigen::Matrix<double, 4, 3> pulley_position_estimate,
                std::vector<double> first_cable_force_magnitude,
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
    forward_kinematic_factor_graph_optimizer(cable_offset, cable_length_collection, cable_forces_collection, p_platform_collection, 
                                             rot_init_platform_collection, delta_rot_platform_collection, pulley_position_estimate,
                                             b_in_w_collection, p_in_w_collection,
                                             first_cable_force_magnitude,
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
    Eigen::Vector4<double> offset = optimization_result.at<gtsam::Vector4>(Symbol('o', 0));

    result->offset = offset;
    result->pulley_result = pulley_optimized;
    result->optimizer_error = optimizer_error; 
}
