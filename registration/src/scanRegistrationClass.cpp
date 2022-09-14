//
// Created by tim on 16.02.21.
//

#include "scanRegistrationClass.h"





Eigen::Matrix4d scanRegistrationClass::sofftRegistration2D(pcl::PointCloud<pcl::PointXYZ> &pointCloudInputData1,
                                                           pcl::PointCloud<pcl::PointXYZ> &pointCloudInputData2,
                                                           double &fitnessX, double &fitnessY, double goodGuessAlpha,
                                                           bool debug) {

//    const pcl::PointCloud<pcl::PointXYZ> pointCloudInputData1New(pointCloudInputData1.makeShared());
//    const pcl::PointCloud<pcl::PointXYZ> pointCloudInputData2New(pointCloudInputData2.makeShared());

    return mySofftRegistrationClass.registrationOfTwoPCL2D(pointCloudInputData1, pointCloudInputData2, fitnessX,
                                                           fitnessY, goodGuessAlpha, debug);
}

Eigen::Matrix4d scanRegistrationClass::sofftRegistration2D(pcl::PointCloud<pcl::PointXYZ> &pointCloudInputData1,
                                                           pcl::PointCloud<pcl::PointXYZ> &pointCloudInputData2,
                                                           double &fitnessX, double &fitnessY,
                                                           Eigen::Matrix4d initialGuess, bool useInitialGuess,
                                                           bool debug) {

//    const pcl::PointCloud<pcl::PointXYZ> pointCloudInputData1New(pointCloudInputData1.makeShared());
//    const pcl::PointCloud<pcl::PointXYZ> pointCloudInputData2New(pointCloudInputData2.makeShared());

    return mySofftRegistrationClass.registrationOfTwoPCL2D(pointCloudInputData1, pointCloudInputData2, fitnessX,
                                                           fitnessY, initialGuess, useInitialGuess, debug);
}


double scanRegistrationClass::sofftRegistrationVoxel2DRotationOnly(double voxelData1Input[], double voxelData2Input[],
                                                                   double goodGuessAlpha, bool debug) {


    return mySofftRegistrationClass.sofftRegistrationVoxel2DRotationOnly(voxelData1Input, voxelData2Input,
                                                                         goodGuessAlpha, debug);

}

//Eigen::Vector2d
//scanRegistrationClass::sofftRegistrationVoxel2DTranslation(double voxelData1Input[], double voxelData2Input[],
//                                                           double &fitnessX, double &fitnessY, double cellSize,
//                                                           Eigen::Vector3d initialGuess, bool useInitialGuess,
//                                                           bool debug) {
//
//    return mySofftRegistrationClass.sofftRegistrationVoxel2DTranslation(voxelData1Input, voxelData2Input, fitnessX,
//                                                                        fitnessY, cellSize, initialGuess,
//                                                                        useInitialGuess, debug);
//
//}




Eigen::Vector2d scanRegistrationClass::sofftRegistrationVoxel2DTranslation(double voxelData1Input[],
                                                                           double voxelData2Input[],
                                                                           double &fitnessX, double &fitnessY,
                                                                           double cellSize,
                                                                           Eigen::Vector3d initialGuess,
                                                                           bool useInitialGuess,
                                                                           double &heightMaximumPeak, bool debug) {

    this->mySofftRegistrationClass.sofftRegistrationVoxel2DTranslation(voxelData1Input, voxelData2Input, fitnessX,
                                                                       fitnessY, cellSize, initialGuess,
                                                                       useInitialGuess, heightMaximumPeak, debug);

}

Eigen::Matrix4d scanRegistrationClass::registrationOfTwoVoxelsSOFFTFast(double voxelData1Input[],
                                                                        double voxelData2Input[],
                                                                        Eigen::Matrix4d initialGuess,
                                                                        bool useInitialAngle, bool useInitialTranslation,
                                                                        double cellSize,
                                                                        bool useGauss,
                                                                        bool debug){



    //changing voxel 1 and 2 because we want to have the transformation from 1 to 2 and not from 2 to 1(which is the registration)
    return mySofftRegistrationClass.registrationOfTwoVoxelsSOFFTFast(voxelData2Input,
                                                                     voxelData1Input,
                                                                     initialGuess,
                                                                     useInitialAngle, useInitialTranslation,
                                                                     cellSize,
                                                                     useGauss,
                                                                     debug);
}