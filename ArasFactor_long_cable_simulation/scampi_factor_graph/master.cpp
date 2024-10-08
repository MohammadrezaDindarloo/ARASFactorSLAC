#include "src/main.cpp"

int main(int argc, char *argv[])
{  
    std::vector<gtsam::Vector10> calibration_result;
    int size_of_calib_sample = 1;
    for (int interval = 0; interval < size_of_calib_sample; interval++) 
    {            
        int lenght_of_simulation_data = 50;
        std::default_random_engine generator(std::random_device{}());
        std::uniform_real_distribution<double> distribution_x(-10.0, 10.0);
        std::uniform_real_distribution<double> distribution_y(-10.0, 10.0);
        std::uniform_real_distribution<double> distribution_z(-5.0, 5.0);

        std::uniform_real_distribution<double> distribution_offset(-0.0, 0.0);

        // std::uniform_real_distribution<double> pulley_location_distribution(-0.4/sqrt(3.0), 0.4/sqrt(3.0));
        std::normal_distribution<double> pulley_location_distribution(0.0, 10.0/sqrt(3.0)/3.0);

        // robot characteristic
        CableRobotParams robot_params(0.7100703113867337, 333.54);

        Eigen::Vector3d Pulley_a(-125.0, -110.0, 48.0);
        Eigen::Vector3d Pulley_b( 125.0, -110.0, 48.0);
        Eigen::Vector3d Pulley_c( 125.0,  110.0, 48.0);
        Eigen::Vector3d Pulley_d(-125.0,  110.0, 48.0);

        robot_params.setPulleyPoses(Pulley_a, Pulley_b, Pulley_c, Pulley_d);

        Eigen::Matrix<double, 4, 3> pulley_position_estimate;
        pulley_position_estimate.row(0) = (Eigen::Vector3d (Pulley_a[0] + pulley_location_distribution(generator), Pulley_a[1] + pulley_location_distribution(generator), Pulley_a[2] + pulley_location_distribution(generator)));
        pulley_position_estimate.row(1) = (Eigen::Vector3d (Pulley_b[0] + pulley_location_distribution(generator), Pulley_b[1] + pulley_location_distribution(generator), Pulley_b[2] + pulley_location_distribution(generator)));
        pulley_position_estimate.row(2) = (Eigen::Vector3d (Pulley_c[0] + pulley_location_distribution(generator), Pulley_c[1] + pulley_location_distribution(generator), Pulley_c[2] + pulley_location_distribution(generator)));
        pulley_position_estimate.row(3) = (Eigen::Vector3d (Pulley_d[0] + pulley_location_distribution(generator), Pulley_d[1] + pulley_location_distribution(generator), Pulley_d[2] + pulley_location_distribution(generator)));
        
        // pulley_position_estimate.row(0) = (Eigen::Vector3d (Pulley_a[0] + pulley_location_distribution(generator), Pulley_a[1] + pulley_location_distribution(generator), Pulley_a[2] + pulley_location_distribution(generator)));
        // pulley_position_estimate.row(1) = (Eigen::Vector3d (Pulley_b[0] + 0.0, Pulley_b[1] + 0.0, Pulley_b[2] + 0.0));
        // pulley_position_estimate.row(2) = (Eigen::Vector3d (Pulley_c[0] + 0.0, Pulley_c[1] + 0.0, Pulley_c[2] + 0.0));
        // pulley_position_estimate.row(3) = (Eigen::Vector3d (Pulley_d[0] + 0.0, Pulley_d[1] + 0.0, Pulley_d[2] + 0.0));
        
        std::vector<double> cable_offset =  {distribution_offset(generator), distribution_offset(generator), distribution_offset(generator), distribution_offset(generator)};

        Eigen::Vector3d Ee_a(-0.21 , -0.21 , -0.011);  
        Eigen::Vector3d Ee_b(0.21  , -0.21 , -0.011);
        Eigen::Vector3d Ee_c(0.21  ,  0.21 , -0.011);
        Eigen::Vector3d Ee_d(-0.21 ,  0.21 , -0.011);
        robot_params.setEEAnchors(Ee_a, Ee_b, Ee_c, Ee_d);

        Eigen::Vector3d r_to_cog(0, 0, -0.12);
        robot_params.setCog(r_to_cog);

        gtsam::Rot3 rot_init_;
        double pitch = 0.01 * M_PI/180.0;
        double roll = 0.01 * M_PI/180.0;
        double yaw = 0.01 * M_PI/180.0;
        // double pitch = -1.31 * M_PI/180.0;
        // double roll = 0.68 * M_PI/180.0;
        // double yaw = -0.04 * M_PI/180.0;
        rot_init_ = rot_init_.RzRyRx(yaw, roll, pitch);
        Eigen::Matrix3d rot_init = gtsamRot3ToEigenMatrix(rot_init_);
        // Eigen::Matrix3d rot_init;
        // rot_init << 0.99268615,  0.11337417, -0.04147891,
        //             -0.11309773,  0.99354347,  0.00895918, 
        //             0.04222684, -0.00420248,  0.99909921; 

        std::vector<Eigen::Matrix<double, 4, 1>> cable_length_collection;
        std::vector<Eigen::Matrix<double, 3, 1>> p_platform_collection;
        std::vector<Eigen::Matrix<double, 3, 3>> rot_init_platform_collection;
        std::vector<Eigen::Matrix<double, 3, 3>> delta_rot_platform_collection;
        std::vector<Eigen::Matrix<double, 2, 1>> cable_forces_collection;

        // std::cout << "******** Extracting Dataset ********" << std::endl;
        // // Open the CSV file of real dataset and record them in data vector
        // std::ifstream file_position("./dataset/pos_i_cpp_test.csv"); 
        // std::vector<std::vector<double>> real_data_position;
        // if (file_position) {
        //     std::string line;
        //     while (getline(file_position, line)) {
        //         std::stringstream ss(line);
        //         std::vector<double> row;
        //         std::string val;
        //         while (getline(ss, val, ',')) {
        //             row.push_back(stod(val));
        //         }
        //         real_data_position.push_back(row);
        //     }
        // std::cout << "Number of position data: " << real_data_position.size() << std::endl;
        // } else {
        //     std::cout << "Unable to open file." << std::endl;
        // }
        // double lenght_dataset_for_calibration = real_data_position.size();
        // // Rewrite the data in it's object
        // for (size_t i = 0; i < real_data_position.size(); i++)
        // {
        //     p_platform_collection.push_back(Eigen::Vector3d(real_data_position[i][0], real_data_position[i][1], real_data_position[i][2]));
        // }
        
        // std::ifstream file_delta_orientation("./dataset/delta_R_i_cpp_test.csv");    
        // std::vector<std::vector<double>> real_delta_orientation;
        // if (file_delta_orientation) {
        //     std::string line;
        //     while (getline(file_delta_orientation, line)) {
        //         std::stringstream ss(line);
        //         std::vector<double> row;
        //         std::string val;
        //         while (getline(ss, val, ',')) {
        //             row.push_back(stod(val));
        //         }
        //         real_delta_orientation.push_back(row);
        //     }
        // std::cout << "Number of delta orientation: " << real_delta_orientation.size() << std::endl;
        // } else {
        //     std::cout << "Unable to open file." << std::endl;
        // }
        // // Rewrite the delta in it's object
        // for (size_t i = 0; i < real_delta_orientation.size(); i++)
        // {   
        //     Eigen::Matrix<double, 3, 3> ith_delta_rot_init;
        //     ith_delta_rot_init <<   real_delta_orientation[i][0], real_delta_orientation[i][1], real_delta_orientation[i][2],
        //                             real_delta_orientation[i][3], real_delta_orientation[i][4], real_delta_orientation[i][5],
        //                             real_delta_orientation[i][6], real_delta_orientation[i][7], real_delta_orientation[i][8];
        //     delta_rot_platform_collection.push_back(ith_delta_rot_init);
        // }
    
        // std::ifstream file_orientation("./dataset/R_i_cpp_test.csv");
        // std::vector<std::vector<double>> real_orientation;
        // if (file_orientation) {
        //     std::string line;
        //     while (getline(file_orientation, line)) {
        //         std::stringstream ss(line);
        //         std::vector<double> row;
        //         std::string val;
        //         while (getline(ss, val, ',')) {
        //             row.push_back(stod(val));
        //         }
        //         real_orientation.push_back(row);
        //     }
        // std::cout << "Number of orientation data: " << real_orientation.size() << std::endl;
        // } else {
        //     std::cout << "Unable to open file." << std::endl;
        // }
        // // Rewrite the data in it's object
        // for (size_t i = 0; i < real_orientation.size(); i++)
        // {   
        //     Eigen::Matrix<double, 3, 3> ith_rot_init;
        //     ith_rot_init << real_orientation[i][0], real_orientation[i][1], real_orientation[i][2],
        //                     real_orientation[i][3], real_orientation[i][4], real_orientation[i][5],
        //                     real_orientation[i][6], real_orientation[i][7], real_orientation[i][8];
        //     rot_init_platform_collection.push_back(ith_rot_init);
        // }

        // std::ifstream file_lcat("./dataset/lc_meas_cpp_test.csv");
        // std::vector<std::vector<double>> real_data_lcat;
        // if (file_lcat) {
        //     std::string line;
        //     while (getline(file_lcat, line)) {
        //         std::stringstream ss(line);
        //         std::vector<double> row;
        //         std::string val;
        //         while (getline(ss, val, ',')) {
        //             row.push_back(stod(val));
        //         }
        //         real_data_lcat.push_back(row);
        //     }
        // std::cout << "Number of encoder data: " << real_data_lcat.size() << std::endl;
        // } else {
        //     std::cout << "Unable to open file." << std::endl;
        // }
        // // Rewrite the data in it's object
        // for (size_t i = 0; i < real_data_lcat.size(); i++)
        // {   
        //     cable_length_collection.push_back(Eigen::Vector4d(real_data_lcat[i][0], real_data_lcat[i][1],
        //                                                     real_data_lcat[i][2], real_data_lcat[i][3]));
        // }

        // std::ifstream file_forces("./dataset/forces_cpp_test.csv");
        // std::vector<std::vector<double>> real_data_forces;
        // if (file_forces) {
        //     std::string line;
        //     while (getline(file_forces, line)) {
        //         std::stringstream ss(line);
        //         std::vector<double> row;
        //         std::string val;
        //         while (getline(ss, val, ',')) {
        //             row.push_back(stod(val));
        //         }
        //         real_data_forces.push_back(row);IKresults
        //     }
        // std::cout << "Number of force sensor data: " << real_data_forces.size() << std::endl;
        // } else {
        //     std::cout << "Unable to open file." << std::endl;
        // }
        // // Rewrite the data in it's object
        // for (size_t i = 0; i < real_data_forces.size(); i++)
        // {   
        //     cable_forces_collection.push_back(Eigen::Vector2d(real_data_forces[i][0], real_data_forces[i][1]));
        // }

        // start inverse optimization for data generation
        for (size_t i = 0; i < lenght_of_simulation_data; i++)
        {                       
            auto p_platform = Eigen::Vector3d((0.0) + distribution_x(generator), (0.0) + distribution_y(generator), (15.0) + distribution_z(generator));
            // auto p_platform = Eigen::Vector3d(0.0, 0.0, 20.0);
            std::vector<MatrixXd> IKresults = IK_Factor_Graph_Optimization(robot_params, rot_init, p_platform);
            // std::cout << std::endl << "rot_platform: " << std::endl << IKresults[0] << std::endl;
            // gtsam::Rot3 rotplat = EigenMatrixToGtsamRot3(IKresults[0]);
            // std::cout << std::endl << "rot_platform roll: " << std::endl << rotplat.roll() << std::endl;
            // std::cout << std::endl << "rot_platform pitch: " << std::endl << rotplat.pitch() << std::endl;
            // std::cout << std::endl << "rot_platform yaw: " << std::endl << rotplat.yaw() << std::endl;
            // std::cout << std::endl << "p_platform: " << std::endl << p_platform << std::endl;
            // std::cout << std::endl << "l_cat: " << std::endl << IKresults[1] << std::endl;
            // std::cout << std::endl << "cable_forces: " << std::endl << IKresults[2] << std::endl;
            // std::cout << std::endl << "c1: " << std::endl << IKresults[3] << std::endl;
            // std::cout << std::endl << "c2: " << std::endl << IKresults[4] << std::endl;
            // std::cout << std::endl << "b_in_w: " << std::endl << IKresults[5] << std::endl;

            // std::cout << std::endl << "sagging_1: " << std::endl << IKresults[1].col(0)[0] - (IKresults[5].col(0) - Pulley_a).norm() << std::endl;
            // std::cout << std::endl << "sagging_2: " << std::endl << IKresults[1].col(0)[1] - (IKresults[5].col(1) - Pulley_b).norm() << std::endl;
            // std::cout << std::endl << "sagging_3: " << std::endl << IKresults[1].col(0)[2] - (IKresults[5].col(2) - Pulley_c).norm() << std::endl;
            // std::cout << std::endl << "sagging_4: " << std::endl << IKresults[1].col(0)[3] - (IKresults[5].col(3) - Pulley_d).norm() << std::endl;
            
            int force_positive = 0;        
            for(size_t i = 0; i < 4; i++)
            {
                if (IKresults[2].col(i)[0]>0 && IKresults[2].col(i)[1]>0)
                {
                    force_positive +=1;
                }
            }
            if (force_positive == 4)
            {
                p_platform_collection.push_back(p_platform);
                rot_init_platform_collection.push_back(rot_init);
                delta_rot_platform_collection.push_back(rot_init.inverse() * IKresults[0]);
                cable_length_collection.push_back(IKresults[1]);
                cable_forces_collection.push_back(Eigen::Matrix<double, 2, 1>(IKresults[2].col(0)));            
            }
        }

        std::vector<Eigen::Matrix<double, 5, 1>> pulley_perturbation_result;
        std::vector<gtsam::Pose3> Optimized_pose;
        std::vector<gtsam::Pose3> GT_pose;
        // start forward optimization
        std::vector<MatrixXd> FKresults = FK_Factor_Graph_Optimization(robot_params, cable_offset, cable_length_collection, cable_forces_collection, p_platform_collection, rot_init_platform_collection, delta_rot_platform_collection, pulley_position_estimate, &Optimized_pose, &GT_pose);

        std::cout << std::endl << "-----------------Calibration Reults------------------------" << std::endl;
        double error_pulley_estimated_a = (Eigen::Vector3d(pulley_position_estimate.row(0)) - Pulley_a).norm() * 1000;
        double error_pulley_estimated_b = (Eigen::Vector3d(pulley_position_estimate.row(1)) - Pulley_b).norm() * 1000;
        double error_pulley_estimated_c = (Eigen::Vector3d(pulley_position_estimate.row(2)) - Pulley_c).norm() * 1000;
        double error_pulley_estimated_d = (Eigen::Vector3d(pulley_position_estimate.row(3)) - Pulley_d).norm() * 1000;
        double sum_pulley_error_estimated = error_pulley_estimated_a + error_pulley_estimated_b + error_pulley_estimated_c + error_pulley_estimated_d;
        std::cout << "Pulley initial guess in mm: " << sum_pulley_error_estimated << std::endl;

        double error_pulley_optimized_a = (Eigen::Vector3d(FKresults[0].row(0)) - Pulley_a).norm() * 1000;
        double error_pulley_optimized_b = (Eigen::Vector3d(FKresults[0].row(1)) - Pulley_b).norm() * 1000;
        double error_pulley_optimized_c = (Eigen::Vector3d(FKresults[0].row(2)) - Pulley_c).norm() * 1000;
        double error_pulley_optimized_d = (Eigen::Vector3d(FKresults[0].row(3)) - Pulley_d).norm() * 1000;
        double sum_pulley_error_optimized = error_pulley_optimized_a + error_pulley_optimized_b + error_pulley_optimized_c + error_pulley_optimized_d;
        std::cout << "Pulley A calibration in  mm: " << error_pulley_optimized_a << std::endl;   
        std::cout << "Pulley B calibration in  mm: " << error_pulley_optimized_b << std::endl;   
        std::cout << "Pulley C calibration in  mm: " << error_pulley_optimized_c << std::endl;   
        std::cout << "Pulley D calibration in  mm: " << error_pulley_optimized_d << std::endl;   
        std::cout << "sum of pulley error  after  calibration  in  mm: " << sum_pulley_error_optimized << std::endl;

        double error_offset_a = std::abs(double(FKresults[1](0)) - cable_offset[0]) * 1000;
        double error_offset_b = std::abs(double(FKresults[1](1)) - cable_offset[1]) * 1000;
        double error_offset_c = std::abs(double(FKresults[1](2)) - cable_offset[2]) * 1000;
        double error_offset_d = std::abs(double(FKresults[1](3)) - cable_offset[3]) * 1000;
        double sum_offset_error_optimized = error_offset_a + error_offset_b + error_offset_c + error_offset_d;
        std::cout << "Offset A calibration in  mm: " << error_offset_a << std::endl;   
        std::cout << "Offset B calibration in  mm: " << error_offset_b << std::endl;   
        std::cout << "Offset C calibration in  mm: " << error_offset_c << std::endl;   
        std::cout << "Offset D calibration in  mm: " << error_offset_d << std::endl;   
        std::cout << "sum of offset error  after  calibration  in  mm: " << sum_offset_error_optimized << std::endl;
        std::cout << "Interval: " << interval << std::endl;
        std::cout << "-----------------Calibration Reults------------------------" << std::endl << std::endl;

        calibration_result.push_back({error_pulley_optimized_a, error_pulley_optimized_b, error_pulley_optimized_c, error_pulley_optimized_d, sum_pulley_error_optimized,
                                      error_offset_a, error_offset_b, error_offset_c, error_offset_d, sum_offset_error_optimized});

        std::ofstream file_Optimized_pose_slakc("./result/Optimized_pose_slakc.csv"); // Optimized_pose_slakc      Optimized_pose_localization
        for (const auto& calib : Optimized_pose) // Loop through the vector elements
        {
            file_Optimized_pose_slakc << calib.translation().x() << ',' << calib.translation().y() << ',' << calib.translation().z() << ',' << std::endl; // Write each element, separated by commas, and end the line
        }
        file_Optimized_pose_slakc.close(); // Close the file stream 

        std::ofstream file_GT_pose_slakc("./result/GT_pose_slakc.csv"); // GT_pose_slakc     GT_pose_localization
        for (const auto& calib : GT_pose) // Loop through the vector elements
        {
            file_GT_pose_slakc << calib.translation().x() << ',' << calib.translation().y() << ',' << calib.translation().z() << ',' << std::endl; // Write each element, separated by commas, and end the line
        }
        file_GT_pose_slakc.close(); // Close the file stream 

    }
    std::ofstream file("./result/calibration_result.csv"); // Create a file stream object
    for (const auto& calib : calibration_result) // Loop through the vector elements
    {
        file << calib[0] << "," << calib[1] << "," << calib[2] << "," << calib[3] << "," << calib[4] << "," << calib[5] << "," << calib[6] << "," << calib[7] << "," << calib[8] << "," << calib[9] << std::endl; // Write each element, separated by commas, and end the line
    }
    file.close(); // Close the file stream
    return 0;
}