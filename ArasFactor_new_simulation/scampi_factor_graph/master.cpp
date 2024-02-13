#include "src/main.cpp"

int main(int argc, char *argv[])
{   
    int lenght_of_simulation_data = 15;
    std::default_random_engine generator(std::random_device{}());
    std::uniform_real_distribution<double> distribution_x(-0.4, 0.4);
    std::uniform_real_distribution<double> distribution_y(-1.3, 1.3);
    std::uniform_real_distribution<double> distribution_z(-2.0, 2.0);
    std::uniform_real_distribution<double> pulley_location_distribution(-0.288675135, 0.288675135);

    // robot characteristic
    CableRobotParams robot_params(0.1034955, 43.164);

    Eigen::Vector3d Pulley_a(-1.9874742031097412, -8.319656372070312, 8.471846580505371);
    Eigen::Vector3d Pulley_b(2.5193355532036756, -8.388501748709967, 8.469020753679201);
    Eigen::Vector3d Pulley_c(2.717151941069913, 4.774436992746004, 8.364108863330584);
    Eigen::Vector3d Pulley_d(-1.7965602546229, 4.832889384134232, 8.370128714520508);

    robot_params.setPulleyPoses(Pulley_a, Pulley_b, Pulley_c, Pulley_d);

    Eigen::Matrix<double, 4, 3> pulley_position_estimate;
    pulley_position_estimate.row(0) = (Eigen::Vector3d (Pulley_a[0] + pulley_location_distribution(generator), Pulley_a[1] + pulley_location_distribution(generator), Pulley_a[2] + pulley_location_distribution(generator)));
    pulley_position_estimate.row(1) = (Eigen::Vector3d (Pulley_b[0] + pulley_location_distribution(generator), Pulley_b[1] + pulley_location_distribution(generator), Pulley_b[2] + pulley_location_distribution(generator)));
    pulley_position_estimate.row(2) = (Eigen::Vector3d (Pulley_c[0] + pulley_location_distribution(generator), Pulley_c[1] + pulley_location_distribution(generator), Pulley_c[2] + pulley_location_distribution(generator)));
    pulley_position_estimate.row(3) = (Eigen::Vector3d (Pulley_d[0] + pulley_location_distribution(generator), Pulley_d[1] + pulley_location_distribution(generator), Pulley_d[2] + pulley_location_distribution(generator)));

    Eigen::Vector3d Ee_a(-0.21 , -0.21 , -0.011);  
    Eigen::Vector3d Ee_b(0.21  , -0.21 , -0.011);
    Eigen::Vector3d Ee_c(0.21  ,  0.21 , -0.011);
    Eigen::Vector3d Ee_d(-0.21 ,  0.21 , -0.011);
    robot_params.setEEAnchors(Ee_a, Ee_b, Ee_c, Ee_d);

    Eigen::Vector3d r_to_cog(0, 0, -0.12);
    robot_params.setCog(r_to_cog);

    Eigen::Matrix3d rot_init;
    rot_init << 0.99268615,  0.11337417, -0.04147891,
                -0.11309773,  0.99354347,  0.00895918,
                0.04222684, -0.00420248,  0.99909921; 

    std::vector<Eigen::Matrix<double, 4, 1>> cable_length_collection;
    std::vector<Eigen::Matrix<double, 3, 1>> p_platform_collection;
    std::vector<Eigen::Matrix<double, 3, 3>> rot_init_platform_collection;
    std::vector<Eigen::Matrix<double, 3, 3>> delta_rot_platform_collection;
    std::vector<Eigen::Matrix<double, 2, 1>> cable_forces_collection;

    for (size_t i = 0; i < lenght_of_simulation_data; i++)
    {
        p_platform_collection.push_back(Eigen::Vector3d(3.09173747e-01 + distribution_x(generator), -1.83715841e+00 + distribution_y(generator),  2.18367984e+00 + distribution_z(generator)));
        // start inverse optimization
        std::vector<MatrixXd> IKresults = IK_Factor_Graph_Optimization(robot_params, rot_init, p_platform_collection[i]);
        // std::cout << std::endl << "rot_platform: " << std::endl << IKresults[0] << std::endl;
        // std::cout << std::endl << "l_cat: " << std::endl << IKresults[1] << std::endl;
        // std::cout << std::endl << "cable_forces: " << std::endl << IKresults[2] << std::endl;
        // std::cout << std::endl << "c1: " << std::endl << IKresults[3] << std::endl;
        // std::cout << std::endl << "c2: " << std::endl << IKresults[4] << std::endl;
        // std::cout << std::endl << "b_in_w: " << std::endl << IKresults[5] << std::endl;
        // std::cout << std::endl << "sagging_1: " << std::endl << IKresults[1].col(0)[0] - (IKresults[5].col(0) - Pulley_a).norm() << std::endl;
        // std::cout << std::endl << "sagging_2: " << std::endl << IKresults[1].col(0)[1] - (IKresults[5].col(1) - Pulley_b).norm() << std::endl;
        // std::cout << std::endl << "sagging_3: " << std::endl << IKresults[1].col(0)[2] - (IKresults[5].col(2) - Pulley_c).norm() << std::endl;
        // std::cout << std::endl << "sagging_4: " << std::endl << IKresults[1].col(0)[3] - (IKresults[5].col(3) - Pulley_d).norm() << std::endl;

        rot_init_platform_collection.push_back(rot_init);
        delta_rot_platform_collection.push_back(rot_init.inverse() * IKresults[0]);
        cable_length_collection.push_back(IKresults[1]);
        cable_forces_collection.push_back(Eigen::Matrix<double, 2, 1>(IKresults[2].col(0)));
    }

    std::vector<Eigen::Matrix<double, 5, 1>> pulley_perturbation_result;
    std::vector<gtsam::Pose3> Optimized_pose;
    std::vector<gtsam::Pose3> GT_pose;
    // start forward optimization
    std::vector<MatrixXd> FKresults = FK_Factor_Graph_Optimization(robot_params, cable_length_collection, cable_forces_collection, p_platform_collection, rot_init_platform_collection, delta_rot_platform_collection, pulley_position_estimate, &Optimized_pose, &GT_pose);

    std::cout << std::endl << "-----------------Calibration Reults------------------------" << std::endl;
    double error_pulley_estimated_a = (Eigen::Vector3d(pulley_position_estimate.row(0)) - Pulley_a).norm() * 1000;
    double error_pulley_estimated_b = (Eigen::Vector3d(pulley_position_estimate.row(1)) - Pulley_b).norm() * 1000;
    double error_pulley_estimated_c = (Eigen::Vector3d(pulley_position_estimate.row(2)) - Pulley_c).norm() * 1000;
    double error_pulley_estimated_d = (Eigen::Vector3d(pulley_position_estimate.row(3)) - Pulley_d).norm() * 1000;
    double sum_pulley_error_estimated = error_pulley_estimated_a + error_pulley_estimated_b + error_pulley_estimated_c + error_pulley_estimated_d;
    std::cout << std::endl << "sum of pulley error in initial estimation in mm: " << sum_pulley_error_estimated << std::endl;

    double error_pulley_optimized_a = (Eigen::Vector3d(FKresults[0].row(0)) - Pulley_a).norm() * 1000;
    double error_pulley_optimized_b = (Eigen::Vector3d(FKresults[0].row(1)) - Pulley_b).norm() * 1000;
    double error_pulley_optimized_c = (Eigen::Vector3d(FKresults[0].row(2)) - Pulley_c).norm() * 1000;
    double error_pulley_optimized_d = (Eigen::Vector3d(FKresults[0].row(3)) - Pulley_d).norm() * 1000;
    double sum_pulley_error_optimized = error_pulley_optimized_a + error_pulley_optimized_b + error_pulley_optimized_c + error_pulley_optimized_d;
    std::cout << std::endl << "A of pulley error  after  calibration  in  mm: " << error_pulley_optimized_a << std::endl;   
    std::cout << std::endl << "B of pulley error  after  calibration  in  mm: " << error_pulley_optimized_b << std::endl;   
    std::cout << std::endl << "C of pulley error  after  calibration  in  mm: " << error_pulley_optimized_c << std::endl;   
    std::cout << std::endl << "D of pulley error  after  calibration  in  mm: " << error_pulley_optimized_d << std::endl;   
    std::cout << std::endl << "sum of pulley error  after  calibration  in  mm: " << sum_pulley_error_optimized << std::endl;

    // std::ofstream file_Optimized_pose_slakc("./result/Optimized_pose_slakc.csv"); // Optimized_pose_slakc      Optimized_pose_localization
    // for (const auto& calib : Optimized_pose) // Loop through the vector elements
    // {
    //     file_Optimized_pose_slakc << calib.translation().x() << ',' << calib.translation().y() << ',' << calib.translation().z() << ',' << std::endl; // Write each element, separated by commas, and end the line
    // }
    // file_Optimized_pose_slakc.close(); // Close the file stream 

    // std::ofstream file_GT_pose_slakc("./result/GT_pose_slakc.csv"); // GT_pose_slakc     GT_pose_localization
    // for (const auto& calib : GT_pose) // Loop through the vector elements
    // {
    //     file_GT_pose_slakc << calib.translation().x() << ',' << calib.translation().y() << ',' << calib.translation().z() << ',' << std::endl; // Write each element, separated by commas, and end the line
    // }
    // file_GT_pose_slakc.close(); // Close the file stream 

    return 0;
}