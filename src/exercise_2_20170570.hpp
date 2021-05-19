#pragma once

Eigen::Matrix3d Rw0, R00i, R0i1, R11i, R1i2, R22i, R2i3, R33i, R3ie;
Eigen::Matrix4d Hw0, H00i, H0i1, H11i, H1i2, H22i, H2i3, H33i, H3ie, Hwe;
Eigen::Vector3d Pw0, P00i, P0i1, P11i, P1i2, P22i, P2i3, P33i, P3ie, Pwe, P01;
Eigen::Vector3d P0e, P1e, P2e, P3e;
Eigen::Vector3d rpyw0, rpy00i, rpy0i1, rpy11i, rpy1i2, rpy22i, rpy2i3, rpy33i, rpy3ie;
Eigen::Vector3d w00i, w11i, w22i, w33i;
Eigen::Vector3d Vw0, V01, V12, V23, V3e;

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

Eigen::Matrix4d HomogeneousMatrix (const Eigen::Matrix3d& RotationMatrix, Eigen::Vector3d& PositionVector) {
    Eigen::Matrix4d H;

    H << RotationMatrix, PositionVector,
            0,0,0,1;

    return H;
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

void Position_calculate (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {


    Pw0 << gc(0), gc(1), gc(2);
    P00i << 0, 0, 0;
    P0i1 << 0.277, 0.116, 0;
    P11i << 0, 0, 0;
    P1i2 << 0.0635, 0.041, 0;
    P22i << 0, 0, 0;
    P2i3 << 0, 0.109, -0.25;
    P33i << 0, 0, 0;
    P3ie << 0.1, -0.02, -0.32125;

    rpyw0 << 0, 0, 0; /// quaternion into rpy
    rpy0i1 << 0, 0, 0;
    rpy11i << gc(7), 0, 0;
    rpy1i2 << 0, 0, 0;
    rpy22i << 0, gc(8), 0;
    rpy2i3 << 0, 0, 0;
    rpy33i << 0, gc(9), 0;
    rpy3ie << 0, 0, 0;

    Rw0 << rpytoRotation(rpyw0);
    R00i << QuaterniontoRotation(gc.segment(3, 4));
    R0i1 << rpytoRotation(rpy0i1);
    R11i << rpytoRotation(rpy11i);
    R1i2 << rpytoRotation(rpy1i2);
    R22i << rpytoRotation(rpy22i);
    R2i3 << rpytoRotation(rpy2i3);
    R33i << rpytoRotation(rpy33i);
    R3ie << rpytoRotation(rpy3ie);

    Hw0 << HomogeneousMatrix(Rw0, Pw0);
    H00i << HomogeneousMatrix(R00i, P00i);
    H0i1 << HomogeneousMatrix(R0i1, P0i1);
    H11i << HomogeneousMatrix(R11i, P11i);
    H1i2 << HomogeneousMatrix(R1i2, P1i2);
    H22i << HomogeneousMatrix(R22i, P22i);
    H2i3 << HomogeneousMatrix(R2i3, P2i3);
    H33i << HomogeneousMatrix(R33i, P33i);
    H3ie << HomogeneousMatrix(R3ie, P3ie);

    Pwe << (Hw0*H00i*H0i1*H11i*H1i2*H22i*H2i3*H33i*H3ie)(12), (Hw0*H00i*H0i1*H11i*H1i2*H22i*H2i3*H33i*H3ie)(13), (Hw0*H00i*H0i1*H11i*H1i2*H22i*H2i3*H33i*H3ie)(14);
    P0e << (H00i*H0i1*H11i*H1i2*H22i*H2i3*H33i*H3ie)(12), (H00i*H0i1*H11i*H1i2*H22i*H2i3*H33i*H3ie)(13), (H00i*H0i1*H11i*H1i2*H22i*H2i3*H33i*H3ie)(14);
    P1e << (H11i*H1i2*H22i*H2i3*H33i*H3ie)(12), (H11i*H1i2*H22i*H2i3*H33i*H3ie)(13), (H11i*H1i2*H22i*H2i3*H33i*H3ie)(14);
    P2e << (H22i*H2i3*H33i*H3ie)(12), (H22i*H2i3*H33i*H3ie)(13), (H22i*H2i3*H33i*H3ie)(14);
    P3e << (H33i*H3ie)(12), (H33i*H3ie)(13), (H33i*H3ie)(14);


    w00i << gv(3), gv(4), gv(5);
    w11i << gv(6), 0, 0;
    w22i << 0, gv(7), 0;
    w33i << 0, gv(8), 0;
}


Eigen::Vector3d CrossProduct (const Eigen::Vector3d& A, Eigen::Matrix<double, 3, 1> B) {
    Eigen::Vector3d Cross;

    Cross << (A(1)*B(2)-A(2)*B(1)), (A(2)*B(0)-A(0)*B(2)), (A(0)*B(1)-A(1)*B(0));

    return Cross;
}

/// do not change the name of the method
inline Eigen::Vector3d getFootLinearVelocity (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
    Eigen::Vector3d LinearVelocity;

    Position_calculate(gc, gv);

    Vw0 << gv(0), gv(1), gv(2);
    V01 = CrossProduct(Rw0*w00i, Rw0*P0e);
    V12 = CrossProduct(Rw0*R00i*R0i1*w11i, Rw0*R00i*R0i1*P1e);
    V23 = CrossProduct(Rw0*R00i*R0i1*R11i*R1i2*w22i, Rw0*R00i*R0i1*R11i*R1i2*P2e);
    V3e = CrossProduct(Rw0*R00i*R0i1*R11i*R1i2*R22i*R2i3*w33i, Rw0*R00i*R0i1*R11i*R1i2*R22i*R2i3*P3e);

    LinearVelocity = Vw0 + V01 + V12 + V23 + V3e;


    return LinearVelocity; /// replace this
}

/// do not change the name of the method
inline Eigen::Vector3d getFootAngularVelocity (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
    Eigen::Vector3d AngularVelocity;
    Position_calculate(gc, gv);

    AngularVelocity = Rw0*(w00i)+ Rw0*R00i*R0i1*(w11i) + Rw0*R00i*R0i1*R11i*R1i2*w22i + Rw0*R00i*R0i1*R11i*R1i2*R22i*R2i3*w33i;

  return AngularVelocity; /// replace this
}