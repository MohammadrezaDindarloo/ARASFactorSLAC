#include "src/main.cpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <filesystem>
#include <random>
#include <bits/stdc++.h>
#include <chrono>





void saveData(string fileName, std::vector<MatrixXd> input_mat)
{
    // Open the CSV file for writing
    std::ofstream csvFile(fileName);
    if (!csvFile.is_open()) {
        std::cerr << "Error opening CSV file." << std::endl;
    }
    // Iterate through the list of matrices
    for (const auto& matrix : input_mat) {
        // Iterate through matrix rows
        for (int i = 0; i < matrix.rows(); ++i) {
            // Iterate through matrix columns
            for (int j = 0; j < matrix.cols(); ++j) {
                // Write matrix element to CSV file
                csvFile << matrix(i, j);

                // Separate values with a comma, except for the last column
                if (j < matrix.cols() - 1) {
                    csvFile << ",";
                }
            }
            
            // Add line break at the end of the row
            csvFile << std::endl;
        }

        // Add an empty line between matrices
        csvFile << std::endl;
    }

    // Close the CSV file
    csvFile.close();
}

void deleteDirectoryContents(const std::filesystem::path& dir)
{
    for (const auto& entry : std::filesystem::directory_iterator(dir)) 
        std::filesystem::remove_all(entry.path());
}

int main(int argc, char *argv[])
{   
    // robot characteristic
    CableRobotParams robot_params(0.1034955, 43.164);

    Eigen::Vector3d Pulley_a(-1.9874742 , -8.31965637,  8.47184658);
    Eigen::Vector3d Pulley_b(2.52022147, -8.38887501,  8.46931362);
    Eigen::Vector3d Pulley_c(2.71799795, 4.77520639, 8.36416322);
    Eigen::Vector3d Pulley_d(-1.79662371,  4.83333111,  8.37001991);
    robot_params.setPulleyPoses(Pulley_a, Pulley_b, Pulley_c, Pulley_d);

    Eigen::Vector3d Ee_a(-0.21 , -0.21 , -0.011);  
    Eigen::Vector3d Ee_b(0.21  , -0.21 , -0.011);
    Eigen::Vector3d Ee_c(0.21  ,  0.21 , -0.011);
    Eigen::Vector3d Ee_d(-0.21 ,  0.21 , -0.011);
    robot_params.setEEAnchors(Ee_a, Ee_b, Ee_c, Ee_d);

    Eigen::Vector3d r_to_cog(0, 0, -0.12);
    robot_params.setCog(r_to_cog);

    Eigen::Vector3d p_platform(3.09173747e-01, -1.83715841e+00,  2.18367984e+00);
    Eigen::Matrix3d rot_init;
    rot_init << 0.99268615,  0.11337417, -0.04147891,
               -0.11309773,  0.99354347,  0.00895918,
                0.04222684, -0.00420248,  0.99909921;

    // read the csv for R_i_cpp_test
    std::ifstream file1("raw_python_data_fedding_cpp_optimization/R_i_cpp_test.csv");
    std::vector<std::vector<double>> R_i_cpp_test;
    if (file1) {
        std::string line;
        while (getline(file1, line)) {
            std::stringstream ss(line);
            std::vector<double> row;
            std::string val;
            while (getline(ss, val, ',')) {
                row.push_back(stod(val));
            }
            R_i_cpp_test.push_back(row);
        }
    std::cout << "Number of R_i_cpp_test: " << R_i_cpp_test.size() << std::endl;
    } else {
        std::cout << "Unable to open file." << std::endl;
    }

    // read the csv for pos_i_cpp_test
    std::ifstream file2("raw_python_data_fedding_cpp_optimization/pos_i_cpp_test.csv");
    std::vector<std::vector<double>> pos_i_cpp_test;
    if (file2) {
        std::string line;
        while (getline(file2, line)) {
            std::stringstream ss(line);
            std::vector<double> row;
            std::string val;
            while (getline(ss, val, ',')) {
                row.push_back(stod(val));
            }
            pos_i_cpp_test.push_back(row);
        }
    std::cout << "Number of pos_i_cpp_test: " << pos_i_cpp_test.size() << std::endl;
    } else {
        std::cout << "Unable to open file." << std::endl;
    }

    // read the csv for lc_meas_cpp_test
    std::ifstream file3("raw_python_data_fedding_cpp_optimization/lc_meas_cpp_test.csv");
    std::vector<std::vector<double>> lc_meas_cpp_test;
    if (file3) {
        std::string line;
        while (getline(file3, line)) {
            std::stringstream ss(line);
            std::vector<double> row;
            std::string val;
            while (getline(ss, val, ',')) {
                row.push_back(stod(val));
            }
            lc_meas_cpp_test.push_back(row);
        }
    std::cout << "Number of lc_meas_cpp_test: " << lc_meas_cpp_test.size() << std::endl;
    } else {
        std::cout << "Unable to open file." << std::endl;
    }


    std::vector<MatrixXd> IKresults_test_rotplatform_IK;
    std::vector<MatrixXd> IKresults_test_lcat_IK;
    std::vector<MatrixXd> IKresults_cableforces_IK;
    std::vector<MatrixXd> IKresults_test_c1_IK;
    std::vector<MatrixXd> IKresults_test_c2_IK;
    std::vector<MatrixXd> IKresults_test_b_in_w_IK;

    std::vector<MatrixXd> FKresults_test_rotplatform_FK;
    std::vector<MatrixXd> FKresults_test_p_platform_FK;
    std::vector<MatrixXd> FKresults_cableforces_FK;
    std::vector<MatrixXd> FKresults_test_c1_FK;
    std::vector<MatrixXd> FKresults_test_c2_FK;
    std::vector<MatrixXd> FKresults_test_b_in_w_FK;

    Eigen::Matrix3d R_rand, R_fk_init;
    Eigen::Vector3d center_fk_init, r_rand;

    auto start = chrono::high_resolution_clock::now();
    // unsync the I/O of C and C++.
    ios_base::sync_with_stdio(false);

    // run the optimization and save the results
    for (size_t i = 0; i < pos_i_cpp_test.size(); i++)
    {   
        Eigen::Vector3d p_platform;
        p_platform[0] = pos_i_cpp_test[i][0];
        p_platform[1] = pos_i_cpp_test[i][1];
        p_platform[2] = pos_i_cpp_test[i][2];
        Eigen::Matrix3d rot_init;
        rot_init << R_i_cpp_test[i][0],  R_i_cpp_test[i][1], R_i_cpp_test[i][2],
                    R_i_cpp_test[i][3],  R_i_cpp_test[i][4], R_i_cpp_test[i][5],
                    R_i_cpp_test[i][6],  R_i_cpp_test[i][7], R_i_cpp_test[i][8];   
        std::vector<MatrixXd> IKresults = IK_Factor_Graph_Optimization(robot_params, rot_init, p_platform);
        IKresults_test_rotplatform_IK.push_back(IKresults[0]);
        IKresults_test_lcat_IK.push_back(IKresults[1]);
        IKresults_cableforces_IK.push_back(IKresults[2]);
        IKresults_test_c1_IK.push_back(IKresults[3]);
        IKresults_test_c2_IK.push_back(IKresults[4]);
        IKresults_test_b_in_w_IK.push_back(IKresults[5]);
        
        R_rand = Eigen::Matrix3d::Identity() + Eigen::Matrix3d::Random() * (0.5 * M_PI / 180.0);
        R_fk_init = IKresults[0] ; //* R_rand

        r_rand = Eigen::Vector3d::Random() * 0.1;
        center_fk_init = p_platform ; //+ r_rand

        Eigen::Vector4d lc_cat;
        lc_cat[0] = lc_meas_cpp_test[i][0]; 
        lc_cat[1] = lc_meas_cpp_test[i][1]; 
        lc_cat[2] = lc_meas_cpp_test[i][2]; 
        lc_cat[3] = lc_meas_cpp_test[i][3];

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dis(0.95, 1.05);
        double rand_value = dis(gen);
        Eigen::Vector2d fc_1 = IKresults[2].col(0) ; //* (rand_value)
        Eigen::Vector3d pos_init = center_fk_init;
        Eigen::Matrix3d rtation_init = R_fk_init;
        std::vector<MatrixXd> FKresults = FK_Factor_Graph_Optimization(robot_params, lc_cat, fc_1, pos_init, rtation_init);
        FKresults_test_rotplatform_FK.push_back(FKresults[0]);
        FKresults_test_p_platform_FK.push_back(FKresults[1]);
        FKresults_cableforces_FK.push_back(FKresults[2]);
        FKresults_test_c1_FK.push_back(FKresults[3]);
        FKresults_test_c2_FK.push_back(FKresults[4]);
        FKresults_test_b_in_w_FK.push_back(FKresults[5]);
        
    }

    auto end = chrono::high_resolution_clock::now();
    double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;
    cout << "Time taken by program is : " << fixed << time_taken << setprecision(9);
    cout << " sec" << endl;

    // clean the directory to save data
    deleteDirectoryContents("data_from_cpp_to_python");

    // save data for using in jupyter notebook
    saveData("data_from_cpp_to_python/IKresults_test_rotplatform_IK.csv", IKresults_test_rotplatform_IK);
    saveData("data_from_cpp_to_python/IKresults_test_lcat_IK.csv", IKresults_test_lcat_IK);
    saveData("data_from_cpp_to_python/IKresults_cableforces_IK.csv", IKresults_cableforces_IK);
    saveData("data_from_cpp_to_python/IKresults_test_c1_IK.csv", IKresults_test_c1_IK);
    saveData("data_from_cpp_to_python/IKresults_test_c2_IK.csv", IKresults_test_c2_IK);
    saveData("data_from_cpp_to_python/IKresults_test_b_in_w_IK.csv", IKresults_test_b_in_w_IK);

    saveData("data_from_cpp_to_python/FKresults_test_rotplatform_FK.csv", FKresults_test_rotplatform_FK);
    saveData("data_from_cpp_to_python/FKresults_test_p_platform_FK.csv", FKresults_test_p_platform_FK);
    saveData("data_from_cpp_to_python/FKresults_cableforces_FK.csv", FKresults_cableforces_FK);
    saveData("data_from_cpp_to_python/FKresults_test_c1_FK.csv", FKresults_test_c1_FK);
    saveData("data_from_cpp_to_python/FKresults_test_c2_FK.csv", FKresults_test_c2_FK);
    saveData("data_from_cpp_to_python/FKresults_test_b_in_w_FK.csv", FKresults_test_b_in_w_FK);

    return 0;
}
