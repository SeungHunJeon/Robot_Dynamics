#pragma once

/// do not change the name of the method

Eigen::Matrix3d R01, R12, R23, R34, R45, R56, R6e;
Eigen::Matrix3d R02, R03, R04, R05, R06, R0e;
Eigen::Vector3d P1e, P2e, P3e, P4e, P5e;
Eigen::Matrix4d H01, H11i, H1i2, H22i, H2i3, H33i, H3i4, H44i, H4i5, H55i, H5i6, H66i, H6ie, H0e;
Eigen::Matrix4d H1e, H2e, H3e, H4e, H5e, H6e;
Eigen::Vector3d P01, P12, P23, P34, P45, P56, P6e, P0e, zero;
Eigen::Vector3d rpy01, rpy12, rpy23, rpy34, rpy45, rpy56, rpy6e;
Eigen::Matrix<double, 3, 6> Ja, Jv;
Eigen::Matrix<double, 6, 6> J;
Eigen::MatrixXd JJ;
Eigen::Vector3d k;
Eigen::Vector3d position_err;
Eigen::Matrix3d R_err;
Eigen::Matrix3d R_des;
Eigen::Vector4d q_err;
Eigen::Vector3d rot_vec_err;
Eigen::VectorXd err;


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
void calculate (const Eigen::VectorXd& gc) {




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

    R02 << R01*jointangletoRotation(gc(0))*R12;
    R03 << R02*jointangletoRotation(gc(1))*R23;
    R04 << R03*jointangletoRotation(gc(2))*R34;
    R05 << R04*jointangletoRotation(gc(3))*R45;
    R06 << R05*jointangletoRotation(gc(4))*R56;

    H6e << H66i*H6ie;
    H5e << H55i*H55i*H6e;
    H4e << H44i*H4i5*H5e;
    H3e << H33i*H3i4*H4e;
    H2e << H22i*H2i3*H3e;
    H1e << H11i*H1i2*H2e;

    P1e << H1e(12), H1e(13), H1e(14);
    P2e << H2e(12), H2e(13), H2e(14);
    P3e << H3e(12), H3e(13), H3e(14);
    P4e << H4e(12), H4e(13), H4e(14);
    P5e << H5e(12), H5e(13), H5e(14);
    P6e << H6e(12), H6e(13), H6e(14);
    H0e << H01*H11i*H1i2*H22i*H2i3*H33i*H3i4*H44i*H4i5*H55i*H5i6*H66i*H6ie;

    R0e << H0e(0), H0e(4), H0e(8),
    H0e(1), H0e(5), H0e(9),
    H0e(2), H0e(6), H0e(10);

    P0e << H0e(12), H0e(13), H0e(14);
    k << 0, 0, 1;

}

Eigen::Vector3d CrossProduct (const Eigen::Vector3d& A, Eigen::Matrix<double, 3, 1> B) {
    Eigen::Vector3d Cross;

    Cross << (A(1)*B(2)-A(2)*B(1)), (A(2)*B(0)-A(0)*B(2)), (A(0)*B(1)-A(1)*B(0));

    return Cross;
}

Eigen::Matrix3d QuaterniontoRotation (const Eigen::Vector4d& quaternion) {
    Eigen::Matrix3d R;
    Eigen::Vector4d q;
    q = quaternion;
    R << (1-2*q(2)*q(2)-2*q(3)*q(3)), 2*(q(1)*q(2)-q(0)*q(3)), 2*(q(1)*q(3)+q(0)*q(2)),
            2*(q(1)*q(2)+q(0)*q(3)), 1-2*q(1)*q(1)-2*q(3)*q(3), 2*(q(2)*q(3)-q(0)*q(1)),
            2*(q(1)*q(3)-q(0)*q(2)), 2*(q(2)*q(3)+q(0)*q(1)), 1-2*q(1)*q(1)-2*q(2)*q(2);
    return R;
}

Eigen::Vector4d RotationtoQuaternion (const Eigen::Matrix3d& R) {
    double m00, m01, m02, m10, m11, m12, m20, m21, m22;

    double qw, qx, qy, qz;
    Eigen::Vector4d q;

    m00 = R(0);
    m10 = R(1);
    m20 = R(2);
    m01 = R(3);
    m11 = R(4);
    m21 = R(5);
    m02 = R(6);
    m12 = R(7);
    m22 = R(8);
    double tr = m00 + m11 + m22;
    if (tr > 0) {
        float S = sqrt(tr+1.0) * 2; // S=4*qw
        qw = 0.25 * S;
        qx = (m21 - m12) / S;
        qy = (m02 - m20) / S;
        qz = (m10 - m01) / S;
    } else if ((m00 > m11)&(m00 > m22)) {
        float S = sqrt(1.0 + m00 - m11 - m22) * 2; // S=4*qx
        qw = (m21 - m12) / S;
        qx = 0.25 * S;
        qy = (m01 + m10) / S;
        qz = (m02 + m20) / S;
    } else if (m11 > m22) {
        float S = sqrt(1.0 + m11 - m00 - m22) * 2; // S=4*qy
        qw = (m02 - m20) / S;
        qx = (m01 + m10) / S;
        qy = 0.25 * S;
        qz = (m12 + m21) / S;
    } else {
        float S = sqrt(1.0 + m22 - m00 - m11) * 2; // S=4*qz
        qw = (m10 - m01) / S;
        qx = (m02 + m20) / S;
        qy = (m12 + m21) / S;
        qz = 0.25 * S;
    }
    q << qw, qx, qy, qz;

    return q;
}

void getJacobian() {
    Ja << R01*k, R02*k, R03*k, R04*k, R05*k, R06*k;
    Jv << CrossProduct(R01*k, R01*P1e), CrossProduct(R02*k, R01*R12*P2e), CrossProduct(R03*k, R01*R12*R23*P3e), CrossProduct(R04*k, R01*R12*R23*R34*P4e), CrossProduct(R05*k, R01*R12*R23*R34*R45*P5e), CrossProduct(R06*k, R01*R12*R23*R34*R45*R56*P6e);

    J << Jv,
    Ja;
}

Eigen::Vector3d quattorotVec(const Eigen::Vector4d& quat) {
    double theta;
    Eigen::Vector3d rotaxis;
    theta = 2*std::acos(quat(0));
    if(theta < 1e-5)
    {
        return Eigen::Vector3d::Zero();
    }

    rotaxis << quat(1)/sin(theta/2), quat(2)/sin(theta/2), quat(3)/sin(theta/2);


    return theta*rotaxis;
}

void getError(const Eigen::Vector3d& pos, const Eigen::Vector4d& quat) {
    err.setZero(6);
    position_err << (pos - P0e);
    R_des << QuaterniontoRotation(quat);
//    std::cout << "R_des : " << R_des << std::endl;
    R_err << R0e.inverse()*R_des;
    q_err << RotationtoQuaternion(R_err);
//    std::cout << "q_err : " << q_err << std::endl;
    rot_vec_err << R0e*quattorotVec(q_err);
//    std::cout << "position_err : " << position_err << std::endl;
//    std::cout << "rot_vec_err : " << rot_vec_err << std::endl;
    err << position_err, rot_vec_err;
}

inline Eigen::VectorXd getVelocityCommand (const Eigen::VectorXd& gc, const Eigen::Vector3d& pos, const Eigen::Vector4d& quat) {
    Eigen::VectorXd Velocity_Command;
    Velocity_Command.setZero(6);

    calculate(gc);
    getJacobian();
    getError(pos, quat);
    JJ = J.completeOrthogonalDecomposition().pseudoInverse(); /// compute pseudoinverse matrix

    Velocity_Command << 100*JJ*err;

//    std::cout<< "Velocity_Command" << Velocity_Command << std::endl;


//    return Eigen::VectorXd::Zero(6);
    return Velocity_Command;
}
