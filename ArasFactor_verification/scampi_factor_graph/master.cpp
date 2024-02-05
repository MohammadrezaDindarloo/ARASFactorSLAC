#include "src/main.cpp"

int main(int argc, char *argv[])
{   
    std::vector<gtsam::Vector5> calibration_result;
    int size_of_calib_sample = 1;
    for (int interval = 0; interval < size_of_calib_sample; interval++)
    {
        double lenght_dataset_for_calibration = 1;
        // create a random engine with a seed
        std::default_random_engine generator(std::random_device{}());
        // create a uniform distribution between 0 and 2
        std::uniform_real_distribution<double> distribution_x(-0.5, 0.5);
        std::uniform_real_distribution<double> distribution_y(-1.0, 1.0);
        std::uniform_real_distribution<double> distribution_z(-1.0, 1.0);

        // robot characteristic
        CableRobotParams robot_params(0.1034955, 43.164);

        Eigen::Vector3d Pulley_a(-1.9874742031097412, -8.319656372070312, 8.471846580505371);
        Eigen::Vector3d Pulley_b(2.5193355532036756, -8.388501748709967, 8.469020753679201);
        Eigen::Vector3d Pulley_c(2.717151941069913, 4.774436992746004, 8.364108863330584);
        Eigen::Vector3d Pulley_d(-1.7965602546229, 4.832889384134232, 8.370128714520508);
        robot_params.setPulleyPoses(Pulley_a, Pulley_b, Pulley_c, Pulley_d);

        std::vector<Eigen::Matrix<double, 3, 1>> pulley_position_estimate;
        pulley_position_estimate.push_back(Pulley_a);
        pulley_position_estimate.push_back(Pulley_b);
        pulley_position_estimate.push_back(Pulley_c);
        pulley_position_estimate.push_back(Pulley_d);

        Eigen::Vector3d Ee_a(-0.21 , -0.21 , -0.011);  
        Eigen::Vector3d Ee_b(0.21  , -0.21 , -0.011);
        Eigen::Vector3d Ee_c(0.21  ,  0.21 , -0.011);
        Eigen::Vector3d Ee_d(-0.21 ,  0.21 , -0.011);
        robot_params.setEEAnchors(Ee_a, Ee_b, Ee_c, Ee_d);

        Eigen::Vector3d r_to_cog(0, 0, -0.12);
        robot_params.setCog(r_to_cog);

        Eigen::Matrix3d rot_init;
        rot_init << 0.99268615,  0.11337417, -0.04147891,
                    -0.11309773,  0.9354347,  0.00895918,
                    0.04222684, -0.00420248,  0.99909921; 

        std::cout << std::endl << "-------------------inverse result--------------------------" << std::endl;
        std::vector<Eigen::Matrix<double, 4, 1>> cable_length_collection;
        std::vector<Eigen::Matrix<double, 3, 1>> p_platform_collection;
        std::vector<Eigen::Matrix<double, 3, 3>> rot_platform_collection;
        std::vector<Eigen::Matrix<double, 3, 3>> delta_rot_platform_collection;
        std::vector<Eigen::Matrix<double, 2, 1>> cable_forces_collection;
        for (size_t i = 0; i < lenght_dataset_for_calibration; i++)
        {   
            double random_number_x = distribution_x(generator);
            double random_number_y = distribution_y(generator);
            double random_number_z = distribution_z(generator);
            p_platform_collection.push_back(Eigen::Vector3d((0.35)+random_number_x, (-1.8)+random_number_y, (3.5)+random_number_z));
            // start inverse optimization
            std::vector<MatrixXd> IKresults = IK_Factor_Graph_Optimization(robot_params, rot_init, p_platform_collection[i]);
            Eigen::Matrix<double, 3, 3> delta_rot = rot_init.inverse() * IKresults[0];
            // the result of inverse optimizationcd 
            rot_platform_collection.push_back(rot_init);
            delta_rot_platform_collection.push_back(delta_rot);
            cable_length_collection.push_back(IKresults[1]);
            cable_forces_collection.push_back(Eigen::Matrix<double, 2, 1>(IKresults[2].col(0)));
        }

        // start forward optimization
        std::vector<MatrixXd> FKresults = FK_Factor_Graph_Optimization(robot_params, pulley_position_estimate, cable_forces_collection, p_platform_collection, rot_platform_collection, cable_length_collection, delta_rot_platform_collection);
        // the result of forward optimization
        std::cout << std::endl << "-------------------forward result--------------------------" << std::endl;
        std::cout << std::endl << "pulley_position: " << std::endl << FKresults[0] << std::endl;
        // calibration result
        double sum_of_error = (Pulley_a - FKresults[0].col(0)).norm() + (Pulley_b - FKresults[0].col(1)).norm() +
                              (Pulley_c - FKresults[0].col(2)).norm() + (Pulley_d - FKresults[0].col(3)).norm();
        double sum_of_error_pulley_1 = (Pulley_a - FKresults[0].col(0)).norm();
        double sum_of_error_pulley_2 = (Pulley_b - FKresults[0].col(1)).norm();
        double sum_of_error_pulley_3 = (Pulley_c - FKresults[0].col(2)).norm();
        double sum_of_error_pulley_4 = (Pulley_d - FKresults[0].col(3)).norm();

        std::cout << std::endl << "sum  of  calibration error in mm: " << std::endl << sum_of_error * 1000 << std::endl;
        std::cout << std::endl << "pulley 1 calibration error in mm: " << std::endl << sum_of_error_pulley_1 * 1000 << std::endl;
        std::cout << std::endl << "pulley 2 calibration error in mm: " << std::endl << sum_of_error_pulley_2 * 1000 << std::endl;
        std::cout << std::endl << "pulley 3 calibration error in mm: " << std::endl << sum_of_error_pulley_3 * 1000 << std::endl;
        std::cout << std::endl << "pulley 4 calibration error in mm: " << std::endl << sum_of_error_pulley_4 * 1000 << std::endl;
        std::cout << std::endl << "Interval: " << std::endl << interval << std::endl;
        
        calibration_result.push_back({sum_of_error_pulley_1 * 1000, sum_of_error_pulley_2 * 1000, sum_of_error_pulley_3 * 1000, sum_of_error_pulley_4 * 1000, sum_of_error * 1000});

    }

    std::ofstream file("./dataset/simuation_enc_uwb_result.csv"); // Create a file stream object
    for (const auto& calib : calibration_result) // Loop through the vector elements
    {
        file << calib[0] << "," << calib[1] << "," << calib[2] << "," << calib[3] << "," << calib[4] << std::endl; // Write each element, separated by commas, and end the line
    }
    file.close(); // Close the file stream 

    return 0;
}
