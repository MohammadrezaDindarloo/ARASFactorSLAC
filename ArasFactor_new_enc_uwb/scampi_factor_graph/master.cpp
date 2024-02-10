#include "src/main.cpp"

int main(int argc, char *argv[])
{   
    int lenght_of_simulation_data = 60;
    std::default_random_engine generator(std::random_device{}());
    std::uniform_real_distribution<double> distribution_x(-0.5, 0.5);
    std::uniform_real_distribution<double> distribution_y(-1.5, 1.5);
    std::uniform_real_distribution<double> distribution_z(-4.0, 3.0);
    std::uniform_real_distribution<double> pulley_location_distribution(-0.0, 0.0);

    // robot characteristic
    CableRobotParams robot_params(0.1034955, 43.164);

    Eigen::Vector3d Pulley_a(-1.9874742031097412, -8.319656372070312, 8.471846580505371);
    Eigen::Vector3d Pulley_b(2.5193355532036756, -8.388501748709967, 8.469020753679201);
    Eigen::Vector3d Pulley_c(2.717151941069913, 4.774436992746004, 8.364108863330584);
    Eigen::Vector3d Pulley_d(-1.7965602546229, 4.832889384134232, 8.370128714520508);

    robot_params.setPulleyPoses(Pulley_a, Pulley_b, Pulley_c, Pulley_d);

    Eigen::Matrix<double, 4, 3> pulley_position_estimate;
    // pulley_position_estimate.row(0) = (Eigen::Vector3d (Pulley_a[0] + pulley_location_distribution(generator), Pulley_a[1] + pulley_location_distribution(generator), Pulley_a[2] + pulley_location_distribution(generator)));
    // pulley_position_estimate.row(1) = (Eigen::Vector3d (Pulley_b[0] + pulley_location_distribution(generator), Pulley_b[1] + pulley_location_distribution(generator), Pulley_b[2] + pulley_location_distribution(generator)));
    // pulley_position_estimate.row(2) = (Eigen::Vector3d (Pulley_c[0] + pulley_location_distribution(generator), Pulley_c[1] + pulley_location_distribution(generator), Pulley_c[2] + pulley_location_distribution(generator)));
    // pulley_position_estimate.row(3) = (Eigen::Vector3d (Pulley_d[0] + pulley_location_distribution(generator), Pulley_d[1] + pulley_location_distribution(generator), Pulley_d[2] + pulley_location_distribution(generator)));

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
    // start forward optimization
    for (int i = -15; i < 15; i++) // x
    {   
        for (int j = -15; j < 15; j++) // y
        {   
            double disturbance_x = double(i)/100.0;
            double disturbance_y = double(j)/100.0;
            pulley_position_estimate.row(0) = (Eigen::Vector3d (Pulley_a[0]+ disturbance_x, Pulley_a[1]+ disturbance_y, Pulley_a[2]));
            pulley_position_estimate.row(1) = (Eigen::Vector3d (Pulley_b[0], Pulley_b[1], Pulley_b[2]));
            pulley_position_estimate.row(2) = (Eigen::Vector3d (Pulley_c[0], Pulley_c[1], Pulley_c[2]));
            pulley_position_estimate.row(3) = (Eigen::Vector3d (Pulley_d[0], Pulley_d[1], Pulley_d[2]));

            std::vector<MatrixXd> FKresults = FK_Factor_Graph_Optimization(robot_params, cable_length_collection, cable_forces_collection, p_platform_collection, rot_init_platform_collection, delta_rot_platform_collection, pulley_position_estimate);
            // the result of forward optimization
            // std::cout << std::endl << "-------------------forward result--------------------------" << std::endl;
            // std::cout << std::endl << "rot_platform: " << std::endl << FKresults[0] << std::endl;
            // std::cout << std::endl << "p_platform: " << std::endl << FKresults[1] << std::endl;
            // std::cout << std::endl << "cable_forces: " << std::endl << FKresults[2] << std::endl;
            // std::cout << std::endl << "c1: " << std::endl << FKresults[3] << std::endl;
            // std::cout << std::endl << "c2: " << std::endl << FKresults[4] << std::endl;
            // std::cout << std::endl << "b_in_w: " << std::endl << FKresults[5] << std::endl;

            // // std::cout << std::endl << "-----------------Calibration Reults------------------------" << std::endl;
            double error_pulley_estimated_a = (Eigen::Vector3d(pulley_position_estimate.row(0)) - Pulley_a).norm() * 1000;
            double error_pulley_estimated_b = (Eigen::Vector3d(pulley_position_estimate.row(1)) - Pulley_b).norm() * 1000;
            double error_pulley_estimated_c = (Eigen::Vector3d(pulley_position_estimate.row(2)) - Pulley_c).norm() * 1000;
            double error_pulley_estimated_d = (Eigen::Vector3d(pulley_position_estimate.row(3)) - Pulley_d).norm() * 1000;
            double sum_pulley_error_estimated = error_pulley_estimated_a + error_pulley_estimated_b + error_pulley_estimated_c + error_pulley_estimated_d;
            std::cout << std::endl << "sum of pulley error in initial estimation in mm: " << sum_pulley_error_estimated << std::endl;

            double error_pulley_optimized_a = (Eigen::Vector3d(FKresults[6].row(0)) - Pulley_a).norm() * 1000;
            double error_pulley_optimized_b = (Eigen::Vector3d(FKresults[6].row(1)) - Pulley_b).norm() * 1000;
            double error_pulley_optimized_c = (Eigen::Vector3d(FKresults[6].row(2)) - Pulley_c).norm() * 1000;
            double error_pulley_optimized_d = (Eigen::Vector3d(FKresults[6].row(3)) - Pulley_d).norm() * 1000;
            double sum_pulley_error_optimized = error_pulley_optimized_a + error_pulley_optimized_b + error_pulley_optimized_c + error_pulley_optimized_d;
            // std::cout << std::endl << "A of pulley error  after  calibration  in  mm: " << error_pulley_optimized_a << std::endl;   
            // std::cout << std::endl << "B of pulley error  after  calibration  in  mm: " << error_pulley_optimized_b << std::endl;   
            // std::cout << std::endl << "C of pulley error  after  calibration  in  mm: " << error_pulley_optimized_c << std::endl;   
            // std::cout << std::endl << "D of pulley error  after  calibration  in  mm: " << error_pulley_optimized_d << std::endl;   
            std::cout << std::endl << "sum of pulley error  after  calibration  in  mm: " << sum_pulley_error_optimized << std::endl;   

            // std::cout << std::endl << "*************************" << std::endl;
            // std::cout << std::endl << "result: " << std::endl << Vector4d(Pulley_a[0] + disturbance_x, Pulley_a[1] + disturbance_y, Pulley_a[2], FKresults[7](0, 0)) << std::endl;
            // std::cout << std::endl << "*************************" << std::endl;
            pulley_perturbation_result.push_back({Pulley_a[0]+disturbance_x, Pulley_a[1]+disturbance_y, Pulley_a[2], FKresults[7](0, 0), sum_pulley_error_optimized});   
            // pulley_perturbation_result.push_back({FKresults[6].row(1)[0], FKresults[6].row(1)[1], FKresults[6].row(1)[2], FKresults[7](0, 0), sum_pulley_error_optimized});   
            // pulley_perturbation_result.push_back({FKresults[6].row(2)[0], FKresults[6].row(2)[1], FKresults[6].row(2)[2], FKresults[7](0, 0), sum_pulley_error_optimized});   
            // pulley_perturbation_result.push_back({FKresults[6].row(3)[0], FKresults[6].row(3)[1], FKresults[6].row(3)[2], FKresults[7](0, 0), sum_pulley_error_optimized});   
        }
    }
    // Save the vectors to a CSV file
    std::ofstream file_pulley("./dataset/pulley_a_perturbation.csv");
    // Assuming all vectors have the same size (1000 in this case)
    size_t vectorSize_pulley = pulley_perturbation_result.size();
    // Iterate through rows
    for (size_t row = 0; row < vectorSize_pulley; row++) {
        // Write values for p_platform_collection as a single cell
        file_pulley  << pulley_perturbation_result[row][0] << "," << pulley_perturbation_result[row][1] << "," << pulley_perturbation_result[row][2] << "," << pulley_perturbation_result[row][3] << "," << pulley_perturbation_result[row][4];
        file_pulley << "\n";
    }
    file_pulley.close();
    return 0;
}