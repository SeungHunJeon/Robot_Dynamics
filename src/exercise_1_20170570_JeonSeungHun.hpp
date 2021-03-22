#pragma once

#define M_PI 3.14159265358979323846

Eigen::Matrix3d rpytoRotation (const Eigen::Vector3d& rpy) {
    Eigen::Matrix3d RotationMatrix;
    Eigen::Matrix3d roll;
    Eigen::Matrix3d pitch;
    Eigen::Matrix3d yaw;

    roll << 1, 0, 0,
            0, cos(rpy(0)), -sin(rpy(0)),
            0, sin(rpy(0)), cos(rpy(0));

    pitch << cos(rpy(1)), 0, sin(rpy(1)),
            0, 1, 0,
            -sin(rpy(1)), 0, cos(rpy(1));

    yaw << cos(rpy(2)), -sin(rpy(2)), 0,
            sin(rpy(2)), cos(rpy(2)), 0,
            0,0,1;

    RotationMatrix = yaw*pitch*roll;

    return RotationMatrix;
}

Eigen::Matrix3d jointangletoRotation (const float& angle) {
    Eigen::Matrix3d RotationMatrix;

    RotationMatrix << cos(angle), -sin(angle), 0,
            sin(angle), cos(angle), 0,
            0,0,1;

    return RotationMatrix;
}

Eigen::Matrix4d HomogeneousMatrix (const Eigen::Matrix3d& RotationMatrix, Eigen::Vector3d& PositionVector) {
    Eigen::Matrix4d H;

    H << RotationMatrix, PositionVector,
    0,0,0,1;

    return H;
}

/// do not change the name of the method
inline Eigen::Vector3d getEndEffectorPosition (const Eigen::VectorXd& gc) {

    Eigen::Matrix3d R01, R12, R23, R34, R45, R56, R6e;
    Eigen::Matrix4d H01, H11i, H1i2, H22i, H2i3, H33i, H3i4, H44i, H4i5, H55i, H5i6, H66i, H6ie, H0e;
    Eigen::Vector3d P01, P12, P23, P34, P45, P56, P6e, P0e, zero;
    Eigen::Vector3d rpy01, rpy12, rpy23, rpy34, rpy45, rpy56, rpy6e;


    P01 << 0, 0, 0.15675;
    P12 << 0, 0.0016, -0.11875;
    P23 << 0, -0.41, 0;
    P34 << 0, 0.2073, -0.0114;
    P45 << 0, 0, -0.10375;
    P56 << 0, 0.10375, 0;
    P6e << 0, 0, -0.16;
    zero << 0, 0, 0;

    rpy01 << M_PI, 0, M_PI;
    rpy12 << -M_PI/2, 0, -M_PI;
    rpy23 << M_PI, 0, M_PI;
    rpy34 << -M_PI/2, 0, -M_PI;
    rpy45 << M_PI/2, 0, -M_PI;
    rpy56 << -M_PI/2, 0, -M_PI;
    rpy6e << -M_PI, 0, 0;

    R01 << rpytoRotation(rpy01);
    R12 << rpytoRotation(rpy12);
    R23 << rpytoRotation(rpy23);
    R34 << rpytoRotation(rpy34);
    R45 << rpytoRotation(rpy45);
    R56 << rpytoRotation(rpy56);
    R6e << rpytoRotation(rpy6e);

    H01 << HomogeneousMatrix(R01, P01);
    H11i << HomogeneousMatrix((jointangletoRotation(gc(0))), zero);
    H1i2 << HomogeneousMatrix(R12, P12);
    H22i << HomogeneousMatrix((jointangletoRotation(gc(1))), zero);
    H2i3 << HomogeneousMatrix(R23, P23);
    H33i << HomogeneousMatrix((jointangletoRotation(gc(2))), zero);
    H3i4 << HomogeneousMatrix(R34, P34);
    H44i << HomogeneousMatrix((jointangletoRotation(gc(3))), zero);
    H4i5 << HomogeneousMatrix(R45, P45);
    H55i << HomogeneousMatrix((jointangletoRotation(gc(4))), zero);
    H5i6 << HomogeneousMatrix(R56, P56);
    H66i << HomogeneousMatrix((jointangletoRotation(gc(5))), zero);
    H6ie << HomogeneousMatrix(R6e, P6e);



    H0e << H01*H11i*H1i2*H22i*H2i3*H33i*H3i4*H44i*H4i5*H55i*H5i6*H66i*H6ie;

    std::cout << H0e << "\n";

    P0e << H0e(12), H0e(13), H0e(14);

    std::cout << P0e << "\n" ;

    return P0e; /// replace this
}

