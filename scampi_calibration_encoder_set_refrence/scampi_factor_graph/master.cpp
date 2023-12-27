#include "src/main.cpp"

int main(int argc, char *argv[])
{   
    double lenght_dataset_for_calibration = 500;

    // create a random engine with a seed
    std::default_random_engine generator(std::random_device{}());
    // create a uniform distribution between 0 and 2
    std::uniform_real_distribution<double> distribution_x(-1.0, 1.0);
    std::uniform_real_distribution<double> distribution_y(-1.0, 1.0);
    std::uniform_real_distribution<double> distribution_z(-1.0, 1.0);

    // robot characteristic
    CableRobotParams robot_params(0.1034955, 43.164);

    Eigen::Vector3d Pulley_a(-1.9874742 , -8.31965637,  8.47184658);
    Eigen::Vector3d Pulley_b(2.52022147, -8.38887501,  8.46931362);
    Eigen::Vector3d Pulley_c(2.71799795, 4.77520639, 8.36416322);
    Eigen::Vector3d Pulley_d(-1.79662371,  4.83333111,  8.37001991);
    robot_params.setPulleyPoses(Pulley_a, Pulley_b, Pulley_c, Pulley_d);

    std::vector<Eigen::Matrix<double, 3, 1>> pulley_position_estimate;
    pulley_position_estimate.push_back(Eigen::Vector3d (-2.5, -8.0,  8.0));
    pulley_position_estimate.push_back(Eigen::Vector3d ( 3.0, -8.0,  8.0));
    pulley_position_estimate.push_back(Eigen::Vector3d ( 3.0,  5.0,  8.0));
    pulley_position_estimate.push_back(Eigen::Vector3d (-2.5,  5.0,  8.0));

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

    // if true, save this collection of data
    if (false)
    {
        std::cout << std::endl << "-------------------inverse result--------------------------" << std::endl;
        std::vector<Eigen::Matrix<double, 3, 1>> p_platform_collection;
        std::vector<Eigen::Matrix<double, 3, 3>> rot_platform_collection;
        std::vector<Eigen::Matrix<double, 2, 1>> cable_forces_collection;
        for (size_t i = 0; i < lenght_dataset_for_calibration; i++)
        {   
            double random_number_x = distribution_x(generator);
            double random_number_y = distribution_y(generator);
            double random_number_z = distribution_z(generator);
            p_platform_collection.push_back(Eigen::Vector3d((0.2)+random_number_x, (-1.5)+random_number_y, (3.5)+random_number_z));
            // start inverse optimization
            std::vector<MatrixXd> IKresults = IK_Factor_Graph_Optimization(robot_params, rot_init, p_platform_collection[i]);
            // the result of inverse optimizationcd
            rot_platform_collection.push_back(IKresults[0]);
            cable_forces_collection.push_back(Eigen::Matrix<double, 2, 1>(IKresults[2].col(0)));
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

    std::vector<Eigen::Matrix<double, 3, 1>> p_platform_collection_dataset;
    std::vector<Eigen::Matrix<double, 3, 3>> rot_platform_collection_dataset;
    std::vector<Eigen::Matrix<double, 2, 1>> cable_forces_collection_dataset;
    for (size_t i = 0; i < lenght_dataset_for_calibration; i++)
    {
        p_platform_collection_dataset.push_back(Eigen::Matrix<double, 3, 1>{data[i][0], data[i][1], data[i][2]});
        cable_forces_collection_dataset.push_back(Eigen::Matrix<double, 2, 1>{data[i][12], data[i][13]});
        rot_platform_collection_dataset.push_back(Eigen::Matrix<double, 3, 3>{ {data[i][3], data[i][6], data[i][9]}, {data[i][4], data[i][7], data[i][10]}, {data[i][5], data[i][8], data[i][11]} });
    }
    
    // start forward optimization
    std::vector<MatrixXd> FKresults = FK_Factor_Graph_Optimization(robot_params, pulley_position_estimate, cable_forces_collection_dataset, p_platform_collection_dataset, rot_platform_collection_dataset);
    // std::vector<MatrixXd> FKresults = FK_Factor_Graph_Optimization(robot_params, pulley_position_estimate, cable_forces_collection, p_platform_collection, rot_platform_collection);
    // the result of forward optimization
    std::cout << std::endl << "-------------------forward result--------------------------" << std::endl;
    std::cout << std::endl << "pulley_position: " << std::endl << FKresults[0] << std::endl;
    // calibration result
    double error = std::pow((Pulley_a - FKresults[0].col(0)).norm(),2) + std::pow((Pulley_b - FKresults[0].col(1)).norm(),2) +
                   std::pow((Pulley_c - FKresults[0].col(2)).norm(),2) + std::pow((Pulley_d - FKresults[0].col(3)).norm(),2);
    std::cout << std::endl << "calibration error in mm: " << std::endl << error * 1000 << std::endl;

    return 0;
}
