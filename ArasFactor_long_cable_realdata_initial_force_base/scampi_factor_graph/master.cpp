#include "src/main.cpp"

int main(int argc, char *argv[])
{  
    std::vector<gtsam::Vector10> calibration_result;
    int size_of_calib_sample = 10;
    for (int interval = 0; interval < size_of_calib_sample; interval++) 
    {            
        std::default_random_engine generator(std::random_device{}());
        std::uniform_real_distribution<double> distribution_x(-10.0, 10.0);
        std::uniform_real_distribution<double> distribution_y(-10.0, 10.0);
        std::uniform_real_distribution<double> distribution_z(-5.0, 5.0);

        std::uniform_real_distribution<double> distribution_offset(0.0, 0.0);

        // std::uniform_real_distribution<double> pulley_location_distribution(-0.4/sqrt(3.0), 0.4/sqrt(3.0));
        std::normal_distribution<double> pulley_location_distribution(0.0, 10.0/sqrt(3.0)/3.0);

        // robot characteristic
        CableRobotParams robot_params(0.7100703113867337, 333.54);
        CableRobotParams robot_params_calibration(0.7100703113867337, 333.54);

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
        
        robot_params_calibration.setPulleyPoses(pulley_position_estimate.row(0), pulley_position_estimate.row(1), pulley_position_estimate.row(2), pulley_position_estimate.row(3));

        std::vector<double> cable_offset =  {distribution_offset(generator), distribution_offset(generator), distribution_offset(generator), distribution_offset(generator)};

        Eigen::Vector3d Ee_a(-0.21 , -0.21 , -0.011);  
        Eigen::Vector3d Ee_b(0.21  , -0.21 , -0.011);
        Eigen::Vector3d Ee_c(0.21  ,  0.21 , -0.011);
        Eigen::Vector3d Ee_d(-0.21 ,  0.21 , -0.011);
        robot_params.setEEAnchors(Ee_a, Ee_b, Ee_c, Ee_d);
        robot_params_calibration.setEEAnchors(Ee_a, Ee_b, Ee_c, Ee_d);

        Eigen::Vector3d r_to_cog(0, 0, -0.12);
        robot_params.setCog(r_to_cog);
        robot_params_calibration.setCog(r_to_cog);

        std::vector<Eigen::Matrix<double, 4, 1>> cable_length_collection;
        std::vector<Eigen::Matrix<double, 3, 1>> p_platform_collection;
        std::vector<Eigen::Matrix<double, 3, 3>> rot_init_platform_collection;
        std::vector<Eigen::Matrix<double, 3, 3>> delta_rot_platform_collection;
        std::vector<Eigen::Matrix<double, 2, 1>> cable_forces_collection;
        std::vector<double> first_cable_force_magnitude;

        // create instance of robot params base on initial guess of pulley for calibration presedure
        RobotParameters<double> params_calibration;
        params_calibration.f_g = robot_params_calibration.f_g_;
        params_calibration.g_c = robot_params_calibration.g_c_;

        params_calibration.ef_points.clear();
        params_calibration.pulleys.clear();

        params_calibration.pulleys.push_back(robot_params_calibration.p1_);
        params_calibration.pulleys.push_back(robot_params_calibration.p2_);
        params_calibration.pulleys.push_back(robot_params_calibration.p3_);
        params_calibration.pulleys.push_back(robot_params_calibration.p4_);

        params_calibration.ef_points.push_back(robot_params_calibration.b1_);
        params_calibration.ef_points.push_back(robot_params_calibration.b2_);
        params_calibration.ef_points.push_back(robot_params_calibration.b3_);
        params_calibration.ef_points.push_back(robot_params_calibration.b4_);

        std::cout << "******** Extracting Dataset ********" << std::endl;
        // Open the CSV file of real dataset and record them in data vector
        std::ifstream file_position("./dataset/pos_i_cpp_test.csv"); 
        std::vector<std::vector<double>> real_data_position;
        if (file_position) {
            std::string line;
            while (getline(file_position, line)) {
                std::stringstream ss(line);
                std::vector<double> row;
                std::string val;
                while (getline(ss, val, ',')) {
                    row.push_back(stod(val));
                }
                real_data_position.push_back(row);
            }
        std::cout << "Number of position data: " << real_data_position.size() << std::endl;
        } else {
            std::cout << "Unable to open file." << std::endl;
        }
        double lenght_dataset_for_calibration = real_data_position.size();
        // Rewrite the data in it's object
        for (size_t i = 0; i < real_data_position.size(); i++)
        {
            p_platform_collection.push_back(Eigen::Vector3d(real_data_position[i][0], real_data_position[i][1], real_data_position[i][2]));
        }

        // ******************************
        // // This is used for euler orientation extraction 
        std::ifstream file_orientation("./dataset/R_i_cpp_test.csv");
        std::vector<std::vector<double>> real_orientation;
        if (file_orientation) {
            std::string line;
            while (getline(file_orientation, line)) {
                std::stringstream ss(line);
                std::vector<double> row;
                std::string val;
                while (getline(ss, val, ',')) {
                    row.push_back(stod(val));
                }
                real_orientation.push_back(row);
            }
        std::cout << "Number of orientation data: " << real_orientation.size() << std::endl;
        } else {
            std::cout << "Unable to open file." << std::endl;
        }
        // Rewrite the data in it's object
        for (size_t i = 0; i < real_orientation.size(); i++)
        {   
            gtsam::Rot3 rot_init_;
            double pitch = real_orientation[i][0] * M_PI/180.0;
            double roll = real_orientation[i][1] * M_PI/180.0;
            double yaw = real_orientation[i][2] * M_PI/180.0;
            rot_init_ = rot_init_.Ypr(yaw, pitch, roll);
            Eigen::Matrix3d rot_init = gtsamRot3ToEigenMatrix(rot_init_);
            rot_init_platform_collection.push_back(rot_init);

            gtsam::Rot3 delta_rot_;
            double pitch_deltaRot = 0.01 * M_PI/180.0;
            double roll_deltaRot = 0.01 * M_PI/180.0;
            double yaw_deltaRot = 0.01 * M_PI/180.0;
            Eigen::Matrix3d deltaRot = gtsamRot3ToEigenMatrix(gtsam::Rot3(delta_rot_.Ypr(yaw_deltaRot, pitch_deltaRot, roll_deltaRot)));
            delta_rot_platform_collection.push_back(deltaRot);
        }

        std::ifstream file_lcat("./dataset/lc_meas_cpp_test.csv");
        std::vector<std::vector<double>> real_data_lcat;
        if (file_lcat) {
            std::string line;
            while (getline(file_lcat, line)) {
                std::stringstream ss(line);
                std::vector<double> row;
                std::string val;
                while (getline(ss, val, ',')) {
                    row.push_back(stod(val));
                }
                real_data_lcat.push_back(row);
            }
        std::cout << "Number of encoder data: " << real_data_lcat.size() << std::endl;
        } else {
            std::cout << "Unable to open file." << std::endl;
        }
        // Rewrite the data in it's object
        for (size_t i = 0; i < real_data_lcat.size(); i++)
        {   
            cable_length_collection.push_back(Eigen::Vector4d(real_data_lcat[i][0], real_data_lcat[i][1],
                                                            real_data_lcat[i][2], real_data_lcat[i][3]));
        }

        std::ifstream file_forces("./dataset/forces_cpp_test.csv");
        std::vector<std::vector<double>> real_data_forces;
        if (file_forces) {
            std::string line;
            while (getline(file_forces, line)) {
                std::stringstream ss(line);
                std::vector<double> row;
                std::string val;
                while (getline(ss, val, ',')) {
                    row.push_back(stod(val));
                }
                real_data_forces.push_back(row);
            }
        std::cout << "Number of force sensor data: " << real_data_forces.size() << std::endl;
        } else {
            std::cout << "Unable to open file." << std::endl;
        }
        // Rewrite the data in it's object
        for (size_t i = 0; i < real_data_forces.size(); i++)
        {   
            first_cable_force_magnitude.push_back(Eigen::Vector2d(real_data_forces[i][0], real_data_forces[i][1]).norm());
            double fh0, fv0;
            computeInitCableForcesCalibration<double>(&fh0, &fv0, p_platform_collection[i], rot_init_platform_collection[i], params_calibration);
            cable_forces_collection.push_back(Eigen::Matrix<double, 2, 1>({std::abs(fh0), std::abs(fv0)})); 
        }

        std::vector<Eigen::Matrix<double, 5, 1>> pulley_perturbation_result;
        std::vector<gtsam::Pose3> Optimized_pose;
        std::vector<gtsam::Pose3> GT_pose;
        // start forward optimization
        std::vector<MatrixXd> FKresults = FK_Factor_Graph_Optimization(robot_params, cable_offset, cable_length_collection, cable_forces_collection, p_platform_collection, rot_init_platform_collection, delta_rot_platform_collection, pulley_position_estimate, first_cable_force_magnitude, &Optimized_pose, &GT_pose);

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