#include "src/main.cpp"

int main(int argc, char *argv[])
{   
    double lenght_dataset_for_calibration = 100;

    // create a random engine with a seed
    std::default_random_engine generator(std::random_device{}());
    // create a uniform distribution between two number
    std::uniform_real_distribution<double> distribution_x(-1.0, 1.0);
    std::uniform_real_distribution<double> distribution_y(-1.0, 1.0);
    std::uniform_real_distribution<double> distribution_z(-1.0, 1.0);
    int a = 2;
    std::uniform_real_distribution<double> pich_purteb(-a, a);
    std::uniform_real_distribution<double> roll_purteb(-a, a);
    std::uniform_real_distribution<double> yaw_purteb(-a, a);

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
    std::vector<Eigen::Matrix<double, 3, 1>> p_platform_collection;
    std::vector<Eigen::Matrix<double, 3, 3>> rot_platform_collection;
    std::vector<Eigen::Matrix<double, 4, 1>> cable_forces_collection;
    for (size_t i = 0; i < lenght_dataset_for_calibration; i++)
    {   
        double random_number_x = distribution_x(generator);
        double random_number_y = distribution_y(generator);
        double random_number_z = distribution_z(generator);
        p_platform_collection.push_back(Eigen::Vector3d((0.2)+random_number_x, (-1.5)+random_number_y, (3.5)+random_number_z));
        // start inverse optimization
        std::vector<MatrixXd> IKresults = IK_Factor_Graph_Optimization(robot_params, rot_init, p_platform_collection[i]);
        // the result of inverse optimizationcd
        double random_pich_purteb = pich_purteb(generator);
        double random_roll_purteb = roll_purteb(generator);
        double random_yaw_purteb = yaw_purteb(generator);
        double euler_angles[3] = {random_pich_purteb, random_roll_purteb, random_yaw_purteb};
        double rotation_matrix_[9];
        Eigen::Matrix3d rotation_matrix_noise;
        ceres::EulerAnglesToRotationMatrix(euler_angles, 3, rotation_matrix_);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                rotation_matrix_noise(i ,j) = rotation_matrix_[i * 3 + j];
            }
        }
        rot_platform_collection.push_back(IKresults[0] * rotation_matrix_noise);
        cable_forces_collection.push_back(Eigen::Matrix<double, 4, 1>(IKresults[2].col(0).norm(), IKresults[2].col(1).norm(), IKresults[2].col(2).norm(), IKresults[2].col(3).norm()));
    }

    // start forward optimization
    std::vector<MatrixXd> FKresults = FK_Factor_Graph_Optimization(robot_params, pulley_position_estimate, cable_forces_collection, p_platform_collection, rot_platform_collection);
    // the result of forward optimization
    std::cout << std::endl << "-------------------forward result--------------------------" << std::endl;
    std::cout << std::endl << "pulley_position: " << std::endl << FKresults[0] << std::endl;
    // calibration result
    double error = std::pow((Pulley_a - FKresults[0].col(0)).norm(),2) + std::pow((Pulley_b - FKresults[0].col(1)).norm(),2) +
                   std::pow((Pulley_c - FKresults[0].col(2)).norm(),2) + std::pow((Pulley_d - FKresults[0].col(3)).norm(),2);
    std::cout << std::endl << "calibration error in mm: " << std::endl << error * 1000 << std::endl;

    return 0;
}
