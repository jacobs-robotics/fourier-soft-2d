//
// Created by tim-linux on 26.03.22.
//

//
// Created by jurobotics on 13.09.21.
//



#include "softDescriptorRegistration.h"

int main(int argc, char **argv) {

    bool debug = true;

    std::string current_exec_name = argv[0]; // Name of the current exec program
    std::vector<std::string> all_args;

    if (argc > 3) {
        // arguments are: first_PCL second_PCL output_dir voxelSize
        // size N is based on image
        all_args.assign(argv + 1, argv + argc);
    } else {
        std::cout << "not enough arguments given" << std::endl;
        exit(-1);
    }

    std::string outputDir = all_args[2];

    int dimensionScan = std::stoi(all_args[3]);
    std::cout << "Voxel size: " << dimensionScan << std::endl;


    softDescriptorRegistration scanRegistrationObject(dimensionScan, dimensionScan / 2, dimensionScan / 2,
                                                      dimensionScan / 2 - 1);;

    pcl::PointCloud<pcl::PointXYZ> scan1;
    pcl::PointCloud<pcl::PointXYZ> scan2;


    pcl::io::loadPCDFile(
            all_args[0],
            scan1);
    pcl::io::loadPCDFile(
            all_args[1],
            scan2);

    Eigen::Matrix4d initialGuess = Eigen::Matrix4d::Identity();

    initialGuess << 1, 0, 0, 15,
           0, 1, -0, 0,
            0, 0, 1, -0,
            -0, 0, -0, 1;

    Eigen::Matrix4d estimatedTransformation = scanRegistrationObject.registrationOfTwoPCL2D(scan1, scan2, initialGuess,
                                                                                            false, false,
                                                                                            outputDir, debug);


    cv::Mat trans_mat = (cv::Mat_<double>(2, 3) << estimatedTransformation(0, 0),
            estimatedTransformation(0, 1),
            estimatedTransformation(0, 3),
            estimatedTransformation(1, 0),
            estimatedTransformation(1, 1),
            estimatedTransformation(1, 3));


    std::cout << "Estimated Transformation:" << std::endl;

    std::cout << estimatedTransformation << std::endl;


    return (0);
}
