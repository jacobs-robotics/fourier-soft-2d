//
// Created by jurobotics on 15.09.21.
//

#include "generalHelpfulTools.h"




double generalHelpfulTools::angleDiff(double first, double second) {//first-second
    return atan2(sin(first - second), cos(first - second));
}


Eigen::Matrix4d generalHelpfulTools::interpolationTwo4DTransformations(Eigen::Matrix4d &transformation1,
                                                                       Eigen::Matrix4d &transformation2, double &t) {
    if (t < 0 || t > 1) {
        std::cout << "t value not between 0 and 1: " << t << std::endl;
        exit(-1);
    }
    Eigen::Vector3d translation1 = transformation1.block<3, 1>(0, 3);
    Eigen::Vector3d translation2 = transformation2.block<3, 1>(0, 3);
    Eigen::Quaterniond rot1(transformation1.block<3, 3>(0, 0));
    Eigen::Quaterniond rot2(transformation2.block<3, 3>(0, 0));

//    Eigen::Quaterniond rotTMP = rot1.inverse()*rot2;
//    Eigen::Quaterniond resultingRot = rotTMP.slerp(t, Eigen::Quaterniond(1,0,0,0));
//    Eigen::Vector3d resultingTranslation = (translation2-translation1)*t;
//    Eigen::Matrix4d resultingTransformation = Eigen::Matrix4d::Identity();
//    resultingTransformation.block<3, 3>(0, 0) = resultingRot.toRotationMatrix();
//    resultingTransformation.block<3, 1>(0, 3) = resultingTranslation;

    Eigen::Quaterniond resultingRot = rot1.slerp(t, rot2);
    Eigen::Vector3d resultingTranslation = translation1 * t + translation2 * (1.0 - t);

    Eigen::Matrix4d resultingTransformation = Eigen::Matrix4d::Identity();
    resultingTransformation.block<3, 3>(0, 0) = resultingRot.toRotationMatrix();
    resultingTransformation.block<3, 1>(0, 3) = resultingTranslation;

    return resultingTransformation;


}


Eigen::Matrix4d
generalHelpfulTools::getTransformationMatrix(Eigen::Vector3d &translation, Eigen::Quaterniond &rotation) {
    Eigen::Matrix4d transformation = Eigen::Matrix4d::Identity();
    transformation.block<3, 1>(0, 3) = translation;
    transformation.block<3, 3>(0, 0) = rotation.toRotationMatrix();
    return transformation;
}

double generalHelpfulTools::weighted_mean(const std::vector<double> &data) {
    double mean = 0.0;

    for (int i = 0; i < data.size(); i++) {
        mean += data[i];
    }
    return mean / double(data.size());
}

void generalHelpfulTools::smooth_curve(const std::vector<double> &input, std::vector<double> &smoothedOutput,
                                       int window_half_width) {

    int window_width = 2 * window_half_width + 1;

    int size = input.size();

    std::vector<double> sample(window_width);

    for (int i = 0; i < size; i++) {

        for (int j = 0; j < window_width; j++) {

            int shifted_index = i + j - window_half_width;
            if (shifted_index < 0) shifted_index = 0;
            if (shifted_index > size - 1) shifted_index = size - 1;
            sample[j] = input[shifted_index];

        }

        smoothedOutput.push_back(generalHelpfulTools::weighted_mean(sample));

    }

}




