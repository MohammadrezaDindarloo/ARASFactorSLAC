#include "src/main.cpp"

int main(int argc, char *argv[])
{   
    std::vector<gtsam::Vector5> calibration_result;
    int size_of_calib_sample = 150;
    for (int interval = 0; interval < size_of_calib_sample; interval++)
    {
        // create a random engine with a seed
        std::default_random_engine generator(std::random_device{}());
        // create a uniform distribution between 0 and 2
        std::uniform_real_distribution<double> distribution_x(-0.5, 0.5);
        std::uniform_real_distribution<double> distribution_y(-1.0, 1.0);
        std::uniform_real_distribution<double> distribution_z(-1.0, 1.0);

        std::uniform_real_distribution<double> pulley_location_distribution(-0.5, 0.5);

        // robot characteristic
        CableRobotParams robot_params(0.1034955, 43.164);

        Eigen::Vector3d Pulley_a(-1.9874742031097412, -8.319656372070312, 8.471846580505371);
        Eigen::Vector3d Pulley_b(2.5193355532036756, -8.388501748709967, 8.469020753679201);
        Eigen::Vector3d Pulley_c(2.717151941069913, 4.774436992746004, 8.364108863330584);
        Eigen::Vector3d Pulley_d(-1.7965602546229, 4.832889384134232, 8.370128714520508);
        robot_params.setPulleyPoses(Pulley_a, Pulley_b, Pulley_c, Pulley_d);

        std::vector<Eigen::Matrix<double, 3, 1>> pulley_position_estimate;

        pulley_position_estimate.push_back(Eigen::Vector3d (Pulley_a[0] + pulley_location_distribution(generator), Pulley_a[1] + pulley_location_distribution(generator), Pulley_a[2] + pulley_location_distribution(generator)));
        pulley_position_estimate.push_back(Eigen::Vector3d (Pulley_b[0] + pulley_location_distribution(generator), Pulley_b[1] + pulley_location_distribution(generator), Pulley_b[2] + pulley_location_distribution(generator)));
        pulley_position_estimate.push_back(Eigen::Vector3d (Pulley_c[0] + pulley_location_distribution(generator), Pulley_c[1] + pulley_location_distribution(generator), Pulley_c[2] + pulley_location_distribution(generator)));
        pulley_position_estimate.push_back(Eigen::Vector3d (Pulley_d[0] + pulley_location_distribution(generator), Pulley_d[1] + pulley_location_distribution(generator), Pulley_d[2] + pulley_location_distribution(generator)));

        // pulley_position_estimate.push_back(Pulley_a);
        // pulley_position_estimate.push_back(Pulley_b);
        // pulley_position_estimate.push_back(Pulley_c);
        // pulley_position_estimate.push_back(Pulley_d);

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

        std::vector<Eigen::Matrix<double, 4, 1>> cable_length_collection_dataset;
        std::vector<Eigen::Matrix<double, 3, 1>> p_platform_collection_dataset;
        std::vector<Eigen::Matrix<double, 3, 3>> rot_platform_collection_dataset;
        std::vector<Eigen::Matrix<double, 2, 1>> cable_forces_collection_dataset;
        
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
        std::cout << "Number of real_data_position: " << real_data_position.size() << std::endl;
        } else {
            std::cout << "Unable to open file." << std::endl;
        }
        double lenght_dataset_for_calibration = real_data_position.size();

        // if true, save this collection of data
        if (false)
        {
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
            std::cout << "Number of real_data_position: " << real_data_position.size() << std::endl;
            } else {
                std::cout << "Unable to open file." << std::endl;
            }

            std::ifstream file_orientation("./dataset/R_i_cpp_test.csv");
            std::vector<std::vector<double>> real_data_orientation;
            if (file_orientation) {
                std::string line;
                while (getline(file_orientation, line)) {
                    std::stringstream ss(line);
                    std::vector<double> row;
                    std::string val;
                    while (getline(ss, val, ',')) {
                        row.push_back(stod(val));
                    }
                    real_data_orientation.push_back(row);
                }
            std::cout << "Number of real_data_orientation: " << real_data_orientation.size() << std::endl;
            } else {
                std::cout << "Unable to open file." << std::endl;
            }

            std::cout << std::endl << "-------------------inverse result--------------------------" << std::endl;
            std::vector<Eigen::Matrix<double, 4, 1>> cable_length_collection;
            std::vector<Eigen::Matrix<double, 3, 1>> p_platform_collection;
            std::vector<Eigen::Matrix<double, 3, 3>> rot_platform_collection;
            std::vector<Eigen::Matrix<double, 2, 1>> cable_forces_collection;
            for (size_t i = 0; i < lenght_dataset_for_calibration; i++)
            {   
                double random_number_x = distribution_x(generator);
                double random_number_y = distribution_y(generator);
                double random_number_z = distribution_z(generator);
                p_platform_collection.push_back(Eigen::Vector3d((0.35)+random_number_x, (-1.8)+random_number_y, (3.5)+random_number_z));
                // p_platform_collection.push_back(Eigen::Vector3d(real_data_position[i][0], real_data_position[i][1], real_data_position[i][2]));
                // start inverse optimization
                std::vector<MatrixXd> IKresults = IK_Factor_Graph_Optimization(robot_params, rot_init, p_platform_collection[i]);
                // the result of inverse optimizationcd
                rot_platform_collection.push_back(IKresults[0]);
                cable_length_collection.push_back(IKresults[1]);
                cable_forces_collection.push_back(Eigen::Matrix<double, 2, 1>(IKresults[2].col(0)));
            }

            size_t numColumns = real_data_position[0].size();

            std::vector<double> minValues(numColumns, std::numeric_limits<double>::max());
            std::vector<double> maxValues(numColumns, std::numeric_limits<double>::lowest());

            // Iterate through each column
            for (size_t col = 0; col < numColumns; ++col) {
                for (const auto& row : real_data_position) {
                    // Update min and max values for the current column
                    minValues[col] = std::min(minValues[col], row[col]);
                    maxValues[col] = std::max(maxValues[col], row[col]);
                }
            }
            // Displaying the boundries
            for (size_t col = 0; col < numColumns; ++col) {
                std::cout << "Column " << col + 1 << " - Min: " << minValues[col] << ", Max: " << maxValues[col] << std::endl;
            }

            // Save the vectors to a CSV file
            std::ofstream file("./dataset/ClibrationDataSet.csv");
            // Assuming all vectors have the same size (1000 in this case)
            size_t vectorSize = lenght_dataset_for_calibration;
            // Iterate through rows
            for (size_t row = 0; row < vectorSize; ++row) {
                // Write values for p_platform_collection as a single cell
                file << p_platform_collection[row](0) << "," << p_platform_collection[row](1) << "," << p_platform_collection[row](2);
                file << ",";  // Separate columns with a comma

                // Write values for rot_platform_collection as a single cell
                for (size_t i = 0; i < 9; ++i) {
                    file << rot_platform_collection[row](i);
                    if (i < 8) {
                        file << ",";
                    }
                }
                file << ",";  // Separate columns with a comma
                // Write values for cable_forces_collection as a single cell
                file << cable_forces_collection[row](0) << "," << cable_forces_collection[row](1);
                file << ",";
                file << cable_length_collection[row](0) << "," << cable_length_collection[row](1) << "," << cable_length_collection[row](2) << "," << cable_length_collection[row](3);
                // Move to the next line after each row
                file << "\n";
            }
            file.close();
        }

        // Open the CSV file of data and record them in data vector
        std::ifstream fileinput("./dataset/ClibrationDataSet.csv");
        std::vector<std::vector<double>> data;
        if (fileinput) {
            std::string line;
            while (getline(fileinput, line)) {
                std::stringstream ss(line);
                std::vector<double> row;
                std::string val;
                while (getline(ss, val, ',')) {
                    row.push_back(stod(val));
                }
                data.push_back(row);
            }
        std::cout << "Number of data: " << data.size() << std::endl;
        } else {
            std::cout << "Unable to open file." << std::endl;
        }

        // using real dataset for calibration
        if(true)
        {   
            // Open the CSV file of real dataset and record them in data vector
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
            std::cout << "Number of real_data_forces: " << real_data_forces.size() << std::endl;
            } else {
                std::cout << "Unable to open file." << std::endl;
            }

            // Open the CSV file of real dataset and record them in data vector
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
            std::cout << "Number of real_data_lcat: " << real_data_lcat.size() << std::endl;
            } else {
                std::cout << "Unable to open file." << std::endl;
            }

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
            std::cout << "Number of real_data_position: " << real_data_position.size() << std::endl;
            } else {
                std::cout << "Unable to open file." << std::endl;
            }

            // Open the CSV file of real dataset and record them in data vector
            std::ifstream file_orientation("./dataset/R_i_cpp_test.csv");
            std::vector<std::vector<double>> real_data_orientation;
            if (file_orientation) {
                std::string line;
                while (getline(file_orientation, line)) {
                    std::stringstream ss(line);
                    std::vector<double> row;
                    std::string val;
                    while (getline(ss, val, ',')) {
                        row.push_back(stod(val));
                    }
                    real_data_orientation.push_back(row);
                }
            std::cout << "Number of real_data_orientation: " << real_data_orientation.size() << std::endl;
            } else {
                std::cout << "Unable to open file." << std::endl;
            }

            lenght_dataset_for_calibration = real_data_orientation.size();
            data.clear();
            for (size_t i = 0; i < lenght_dataset_for_calibration; i++)
            {   
                data.push_back({real_data_position[i][0], real_data_position[i][1], real_data_position[i][2],
                                real_data_orientation[i][0], real_data_orientation[i][3], real_data_orientation[i][6],
                                real_data_orientation[i][1], real_data_orientation[i][4], real_data_orientation[i][7],
                                real_data_orientation[i][2], real_data_orientation[i][5], real_data_orientation[i][8],
                                real_data_forces[i][0], real_data_forces[i][1],
                                real_data_lcat[i][0], real_data_lcat[i][1], real_data_lcat[i][2], real_data_lcat[i][3]});
            }
        }

        for (size_t i = 0; i < lenght_dataset_for_calibration; i++)
        {   
            p_platform_collection_dataset.push_back(Eigen::Matrix<double, 3, 1>{data[i][0], data[i][1], data[i][2]});
            cable_forces_collection_dataset.push_back(Eigen::Matrix<double, 2, 1>{data[i][12], data[i][13]});
            rot_platform_collection_dataset.push_back(Eigen::Matrix<double, 3, 3>{ {data[i][3], data[i][6], data[i][9]}, {data[i][4], data[i][7], data[i][10]}, {data[i][5], data[i][8], data[i][11]} });
            cable_length_collection_dataset.push_back(Eigen::Matrix<double, 4, 1>{data[i][14], data[i][15], data[i][16], data[i][17]});
        }

        // start forward optimization
        std::vector<MatrixXd> FKresults = FK_Factor_Graph_Optimization(robot_params, pulley_position_estimate, cable_forces_collection_dataset, p_platform_collection_dataset, rot_platform_collection_dataset, cable_length_collection_dataset);
        // the result of forward optimization
        std::cout << std::endl << "-------------------forward result--------------------------" << std::endl;
        std::cout << std::endl << "pulley_position: " << std::endl << FKresults[0] << std::endl;
        // calibration result
        double sum_of_error = std::pow((Pulley_a - FKresults[0].col(0)).norm(),2) + std::pow((Pulley_b - FKresults[0].col(1)).norm(),2) +
                    std::pow((Pulley_c - FKresults[0].col(2)).norm(),2) + std::pow((Pulley_d - FKresults[0].col(3)).norm(),2);

        double sum_of_error_pulley_1 = std::pow((Pulley_a - FKresults[0].col(0)).norm(),2);

        double sum_of_error_pulley_2 = std::pow((Pulley_b - FKresults[0].col(1)).norm(),2);

        double sum_of_error_pulley_3 = std::pow((Pulley_c - FKresults[0].col(2)).norm(),2);

        double sum_of_error_pulley_4 = std::pow((Pulley_d - FKresults[0].col(3)).norm(),2);

        std::cout << std::endl << "sum  of  calibration error in mm: " << std::endl << sum_of_error * 1000 << std::endl;
        std::cout << std::endl << "pulley 1 calibration error in mm: " << std::endl << sum_of_error_pulley_1 * 1000 << std::endl;
        std::cout << std::endl << "pulley 2 calibration error in mm: " << std::endl << sum_of_error_pulley_2 * 1000 << std::endl;
        std::cout << std::endl << "pulley 3 calibration error in mm: " << std::endl << sum_of_error_pulley_3 * 1000 << std::endl;
        std::cout << std::endl << "pulley 4 calibration error in mm: " << std::endl << sum_of_error_pulley_4 * 1000 << std::endl;
        std::cout << std::endl << "Interval: " << std::endl << interval << std::endl;
          
        calibration_result.push_back({sum_of_error_pulley_1 * 1000, sum_of_error_pulley_2 * 1000, sum_of_error_pulley_3 * 1000, sum_of_error_pulley_4 * 1000, sum_of_error * 1000});

    }
    std::ofstream file("./dataset/realdata_enc_result.csv"); // Create a file stream object
    for (const auto& calib : calibration_result) // Loop through the vector elements
    {
        file << calib[0] << "," << calib[1] << "," << calib[2] << "," << calib[3] << "," << calib[4] << std::endl; // Write each element, separated by commas, and end the line
    }
    file.close(); // Close the file stream

    return 0;
}
