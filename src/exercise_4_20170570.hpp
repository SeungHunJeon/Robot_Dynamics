#pragma once

Eigen::Matrix3d R01, R12, R23, R34, R45, R56, R6e;
Eigen::Matrix3d R01i, R02i, R03i, R04i, R05i, R06i;
Eigen::Matrix3d R02, R03, R04, R05, R06, R0e;
Eigen::Vector3d P1e, P2e, P3e, P4e, P5e;
Eigen::Matrix4d H01, H11i, H1i2, H22i, H2i3, H33i, H3i4, H44i, H4i5, H55i, H5i6, H66i, H6ie, H0e;
Eigen::Matrix4d H02, H03, H04, H05, H06;
Eigen::Matrix4d H1e, H2e, H3e, H4e, H5e, H6e;
Eigen::Vector3d P01, P12, P23, P34, P45, P56, P6e, P0e, zero;
Eigen::Vector3d P02, P03, P04, P05, P06;
Eigen::Vector3d P01com, P02com, P03com, P04com, P05com, P06com, P0ecom;
Eigen::Vector3d P1com, P2com, P3com, P4com, P5com, P6com, Pecom;
Eigen::Vector3d rpy01, rpy12, rpy23, rpy34, rpy45, rpy56, rpy6e;
Eigen::Matrix<double, 3, 6> J1a, J1v, J2a, J2v, J3a, J3v, J4a, J4v, J5a, J5v, J6a, J6v, Jea, Jev;
Eigen::Matrix<double, 3, 6> J1a2, J1v2, J2a2, J2v2, J3a2, J3v2, J4a2, J4v2, J5a2, J5v2, J6a2, J6v2, Jea2, Jev2;
Eigen::Matrix<double, 3, 6> J1a_dot, J1v_dot, J2a_dot, J2v_dot, J3a_dot, J3v_dot, J4a_dot, J4v_dot, J5a_dot, J5v_dot, J6a_dot, J6v_dot, Jea_dot, Jev_dot;
Eigen::Matrix<double, 6, 6> Je;
Eigen::MatrixXd JJ;
Eigen::Vector3d k;
Eigen::Vector3d position_err;
Eigen::Matrix3d R_err;
Eigen::Matrix3d R_des;
Eigen::Vector4d q_err;
Eigen::Vector3d rot_vec_err;
Eigen::VectorXd err;
Eigen::Matrix3d I1, I2, I3, I4, I5, I6, Ie;
double m1, m2, m3, m4, m5, m6, me;

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

Eigen::Matrix3d jointangletoRotation (const double& angle) {
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

    R01i << R01*jointangletoRotation(gc(0));
    R02i << R02*jointangletoRotation(gc(1));
    R03i << R03*jointangletoRotation(gc(2));
    R04i << R04*jointangletoRotation(gc(3));
    R05i << R05*jointangletoRotation(gc(4));
    R06i << R06*jointangletoRotation(gc(5));

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


    H02 << H01*H11i*H1i2;
    H03 << H02*H22i*H2i3;
    H04 << H03*H33i*H3i4;
    H05 << H04*H44i*H4i5;
    H06 << H05*H55i*H5i6;
    H0e << H06*H66i*H6ie;



    R0e << H0e(0), H0e(4), H0e(8),
            H0e(1), H0e(5), H0e(9),
            H0e(2), H0e(6), H0e(10);

    P02 << H02(12), H02(13), H02(14);
    P03 << H03(12), H03(13), H03(14);
    P04 << H04(12), H04(13), H04(14);
    P05 << H05(12), H05(13), H05(14);
    P06 << H06(12), H06(13), H06(14);
    P0e << H0e(12), H0e(13), H0e(14);

    P1com << 0, -0.002, -0.0605;
    P2com << 0, -0.2065, -0.01;
    P3com << 0, 0.081, -0.0086;
    P4com << 0, 0.0028848942, -0.0541932613;
    P5com << 0, 0.0497208855, -0.0028562765;
    P6com << 0, 0, -0.06;
    Pecom << 0, 0, 0;

    P01com << R01i*P1com+P01;
    P02com << R02i*P2com+P02;
    P03com << R03i*P3com+P03;
    P04com << R04i*P4com+P04;
    P05com << R05i*P5com+P05;
    P06com << R06i*P6com+P06;
    P0ecom << R0e*Pecom+P0e;

    k << 0, 0, 1;

    I1 << 0.00152031725204, 0, 0,
            0, 0.00152031725204, 0,
            0, 0, 0.00059816;
    I2 << 0.010502207991, 0, 0,
            0, 0.000792, 0,
            0, 0, 0.010502207991;
    I3 << 0.00142022431908, 0, 0,
            0, 0.000304335, 0,
            0, 0, 0.00142022431908;
    I4 << 0.0004321316048, 0, 0,
            0, 0.0004321316048, 0,
            0, 0, 9.26e-05;
    I5 << 0.0004321316048, 0, 0,
            0, 9.26e-05, 0,
            0, 0, 0.0004321316048;
    I6 << 0.0004403232387, 0, 0,
            0, 0.0004403232387, 0,
            0, 0, 0.0007416;
    Ie << 0.01, 0, 0,
            0, 0.01, 0,
            0, 0, 0.01;

    ///////
    I1 << R01i*I1*R01i.transpose();
    I2 << R02i*I2*R02i.transpose();
    I3 << R03i*I3*R03i.transpose();
    I4 << R04i*I4*R04i.transpose();
    I5 << R05i*I5*R05i.transpose();
    I6 << R06i*I6*R06i.transpose();
    Ie << R0e*Ie*R0e.transpose();
    //////

    m1 = 0.7477;
    m2 = 0.99;
    m3 = 0.6763;
    m4 = 0.463;
    m5 = 0.463;
    m6 = 1.327;
    me = 0.01;
}

inline Eigen::Vector3d getEndEffectorPosition (const Eigen::VectorXd& gc) {

    calculate(gc);

    return P0ecom; /// replace this
}

Eigen::Vector3d CrossProduct (const Eigen::Vector3d& A, Eigen::Matrix<double, 3, 1> B) {
    Eigen::Vector3d Cross;

    Cross << (A(1)*B(2)-A(2)*B(1)), (A(2)*B(0)-A(0)*B(2)), (A(0)*B(1)-A(1)*B(0));

    return Cross;
}


void getJacobian() {
    Eigen::Matrix3d Zero;
    Zero.setZero();

    Jea << R01*k, R02*k, R03*k, R04*k, R05*k, R06*k;
    Jev << CrossProduct(R01*k, (P0ecom-P01)), CrossProduct(R02*k, (P0ecom-P02)), CrossProduct(R03*k, (P0ecom-P03)), CrossProduct(R04*k, (P0ecom-P04)), CrossProduct(R05*k, (P0ecom-P05)), CrossProduct(R06*k, (P0ecom-P06));

    J6a << R01*k, R02*k, R03*k, R04*k, R05*k, R06*k;
    J6v << CrossProduct(R01*k, (P06com-P01)), CrossProduct(R02*k, (P06com-P02)), CrossProduct(R03*k, (P06com-P03)), CrossProduct(R04*k, (P06com-P04)), CrossProduct(R05*k, (P06com-P05)), CrossProduct(R06*k, (P06com-P06));

    J5a << R01*k, R02*k, R03*k, R04*k, R05*k, Zero*k;
    J5v << CrossProduct(R01*k, (P05com-P01)), CrossProduct(R02*k, (P05com-P02)), CrossProduct(R03*k, (P05com-P03)), CrossProduct(R04*k, (P05com-P04)), CrossProduct(R05*k, (P05com-P05)), Zero*k;

    J4a << R01*k, R02*k, R03*k, R04*k, Zero*k, Zero*k;
    J4v << CrossProduct(R01*k, (P04com-P01)), CrossProduct(R02*k, (P04com-P02)), CrossProduct(R03*k, (P04com-P03)), CrossProduct(R04*k, (P04com-P04)), Zero*k, Zero*k;

    J3a << R01*k, R02*k, R03*k, Zero*k, Zero*k, Zero*k;
    J3v << CrossProduct(R01*k, (P03com-P01)), CrossProduct(R02*k, (P03com-P02)), CrossProduct(R03*k, (P03com-P03)), Zero*k, Zero*k, Zero*k;

    J2a << R01*k, R02*k, Zero*k, Zero*k, Zero*k, Zero*k;
    J2v << CrossProduct(R01*k, (P02com-P01)), CrossProduct(R02*k, (P02com-P02)), Zero*k, Zero*k, Zero*k, Zero*k;

    J1a << R01*k, Zero*k, Zero*k, Zero*k, Zero*k, Zero*k;
    J1v << CrossProduct(R01*k, (P01com-P01)), Zero*k, Zero*k, Zero*k, Zero*k, Zero*k;

    //////////////////////////////////////////////////

    Jea2 << R01*k, R02*k, R03*k, R04*k, R05*k, R06*k;
    Jev2 << CrossProduct(R01*k, (P0e-P01)), CrossProduct(R02*k, (P0e-P02)), CrossProduct(R03*k, (P0e-P03)), CrossProduct(R04*k, (P0e-P04)), CrossProduct(R05*k, (P0e-P05)), CrossProduct(R06*k, (P0e-P06));

    J6a2 << R01*k, R02*k, R03*k, R04*k, R05*k, Zero*k;
    J6v2 << CrossProduct(R01*k, (P06-P01)), CrossProduct(R02*k, (P06-P02)), CrossProduct(R03*k, (P06-P03)), CrossProduct(R04*k, (P06-P04)), CrossProduct(R05*k, (P06-P05)), CrossProduct(R06*k, (P06-P06));

    J5a2 << R01*k, R02*k, R03*k, R04*k, Zero*k, Zero*k;
    J5v2 << CrossProduct(R01*k, (P05-P01)), CrossProduct(R02*k, (P05-P02)), CrossProduct(R03*k, (P05-P03)), CrossProduct(R04*k, (P05-P04)), CrossProduct(R05*k, (P05-P05)), Zero*k;

    J4a2 << R01*k, R02*k, R03*k, Zero*k, Zero*k, Zero*k;
    J4v2 << CrossProduct(R01*k, (P04-P01)), CrossProduct(R02*k, (P04-P02)), CrossProduct(R03*k, (P04-P03)), CrossProduct(R04*k, (P04-P04)), Zero*k, Zero*k;

    J3a2 << R01*k, R02*k,  Zero*k, Zero*k, Zero*k, Zero*k;
    J3v2 << CrossProduct(R01*k, (P03-P01)), CrossProduct(R02*k, (P03-P02)), CrossProduct(R03*k, (P03-P03)), Zero*k, Zero*k, Zero*k;

    J2a2 << R01*k, Zero*k, Zero*k, Zero*k, Zero*k, Zero*k;
    J2v2 << CrossProduct(R01*k, (P02-P01)), CrossProduct(R02*k, (P02-P02)), Zero*k, Zero*k, Zero*k, Zero*k;

    J1a2 << Zero*k, Zero*k, Zero*k, Zero*k, Zero*k, Zero*k;
    J1v2 << CrossProduct(R01*k, (P01-P01)), Zero*k, Zero*k, Zero*k, Zero*k, Zero*k;

}

void getJacobian_dot(const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
    Eigen::Matrix3d Zero;
    Zero.setZero();

    Jea_dot << CrossProduct(J1a*gv,R01*k), CrossProduct(J2a*gv, R02*k), CrossProduct(J3a*gv, R03*k), CrossProduct(J4a*gv, R04*k), CrossProduct(J5a*gv, R05*k), CrossProduct(J6a*gv, R06*k);
    Jev_dot << CrossProduct(CrossProduct(J1a*gv, R01*k),(P0ecom-P01)) + CrossProduct(R01*k, Jev*gv-J1v2*gv), CrossProduct(CrossProduct(J2a*gv, R02*k),(P0ecom-P02)) + CrossProduct(R02*k, Jev*gv-J2v2*gv), CrossProduct(CrossProduct(J3a*gv, R03*k),(P0ecom-P03)) + CrossProduct(R03*k, Jev*gv-J3v2*gv), CrossProduct(CrossProduct(J4a*gv, R04*k),(P0ecom-P04)) + CrossProduct(R04*k, Jev*gv-J4v2*gv), CrossProduct(CrossProduct(J5a*gv, R05*k),(P0ecom-P05)) + CrossProduct(R05*k, Jev*gv-J5v2*gv), CrossProduct(CrossProduct(J6a*gv, R06*k),(P0ecom-P06)) + CrossProduct(R06*k, Jev*gv-J6v2*gv);

    J6a_dot << CrossProduct(J1a*gv,R01*k), CrossProduct(J2a*gv, R02*k), CrossProduct(J3a*gv, R03*k), CrossProduct(J4a*gv, R04*k), CrossProduct(J5a*gv, R05*k), CrossProduct(J6a*gv, R06*k);
    J6v_dot << CrossProduct(CrossProduct(J1a*gv, R01*k),(P06com-P01)) + CrossProduct(R01*k, J6v*gv-J1v2*gv), CrossProduct(CrossProduct(J2a*gv, R02*k),(P06com-P02)) + CrossProduct(R02*k, J6v*gv-J2v2*gv), CrossProduct(CrossProduct(J3a*gv, R03*k),(P06com-P03)) + CrossProduct(R03*k, J6v*gv-J3v2*gv), CrossProduct(CrossProduct(J4a*gv, R04*k),(P06com-P04)) + CrossProduct(R04*k, J6v*gv-J4v2*gv), CrossProduct(CrossProduct(J5a*gv, R05*k),(P06com-P05)) + CrossProduct(R05*k, J6v*gv-J5v2*gv), CrossProduct(CrossProduct(J6a*gv, R06*k),(P06com-P06)) + CrossProduct(R06*k, J6v*gv-J6v2*gv);

    J5a_dot << CrossProduct(J1a*gv,R01*k), CrossProduct(J2a*gv, R02*k), CrossProduct(J3a*gv, R03*k), CrossProduct(J4a*gv, R04*k), CrossProduct(J5a*gv, R05*k), Zero*k;
    J5v_dot << CrossProduct(CrossProduct(J1a*gv, R01*k),(P05com-P01)) + CrossProduct(R01*k, J5v*gv-J1v2*gv), CrossProduct(CrossProduct(J2a*gv, R02*k),(P05com-P02)) + CrossProduct(R02*k, J5v*gv-J2v2*gv), CrossProduct(CrossProduct(J3a*gv, R03*k),(P05com-P03)) + CrossProduct(R03*k, J5v*gv-J3v2*gv), CrossProduct(CrossProduct(J4a*gv, R04*k),(P05com-P04)) + CrossProduct(R04*k, J5v*gv-J4v2*gv), CrossProduct(CrossProduct(J5a*gv, R05*k),(P05com-P05)) + CrossProduct(R05*k, J5v*gv-J5v2*gv), Zero*k;

    J4a_dot << CrossProduct(J1a*gv,R01*k), CrossProduct(J2a*gv, R02*k), CrossProduct(J3a*gv, R03*k), CrossProduct(J4a*gv, R04*k), Zero*k, Zero*k;
    J4v_dot << CrossProduct(CrossProduct(J1a*gv, R01*k),(P04com-P01)) + CrossProduct(R01*k, J4v*gv-J1v2*gv), CrossProduct(CrossProduct(J2a*gv, R02*k),(P04com-P02)) + CrossProduct(R02*k, J4v*gv-J2v2*gv), CrossProduct(CrossProduct(J3a*gv, R03*k),(P04com-P03)) + CrossProduct(R03*k, J4v*gv-J3v2*gv), CrossProduct(CrossProduct(J4a*gv, R04*k),(P04com-P04)) + CrossProduct(R04*k, J4v*gv-J4v2*gv), Zero*k, Zero*k;

    J3a_dot << CrossProduct(J1a*gv,R01*k), CrossProduct(J2a*gv, R02*k), CrossProduct(J3a*gv, R03*k), Zero*k, Zero*k, Zero*k;
    J3v_dot << CrossProduct(CrossProduct(J1a*gv, R01*k),(P03com-P01)) + CrossProduct(R01*k, J3v*gv-J1v2*gv), CrossProduct(CrossProduct(J2a*gv, R02*k),(P03com-P02)) + CrossProduct(R02*k, J3v*gv-J2v2*gv), CrossProduct(CrossProduct(J3a*gv, R03*k),(P03com-P03)) + CrossProduct(R03*k, J3v*gv-J3v2*gv), Zero*k, Zero*k, Zero*k;

    J2a_dot << CrossProduct(J1a*gv,R01*k), CrossProduct(J2a*gv, R02*k), Zero*k, Zero*k, Zero*k, Zero*k;
    J2v_dot << CrossProduct(CrossProduct(J1a*gv, R01*k),(P02com-P01)) + CrossProduct(R01*k, J2v*gv-J1v2*gv), CrossProduct(CrossProduct(J2a*gv, R02*k),(P02com-P02)) + CrossProduct(R02*k, J2v*gv-J2v2*gv), Zero*k, Zero*k, Zero*k, Zero*k;

    J1a_dot << CrossProduct(J1a*gv,R01*k), Zero*k, Zero*k, Zero*k, Zero*k, Zero*k;
    J1v_dot << CrossProduct(CrossProduct(J1a*gv, R01*k),(P01com-P01)) + CrossProduct(R01*k, (J1v*gv)-(J1v2*gv)), Zero*k, Zero*k, Zero*k, Zero*k, Zero*k;

}


/// do not change the name of the method
inline Eigen::MatrixXd getMassMatrix (const Eigen::VectorXd& gc) {

    Eigen::MatrixXd M;
    M.resize(6, 6);
    calculate(gc);
    getJacobian();

    M << J1v.transpose()*m1*J1v+J1a.transpose()*I1*J1a + J2v.transpose()*m2*J2v+J2a.transpose()*I2*J2a + J3v.transpose()*m3*J3v+J3a.transpose()*I3*J3a + J4v.transpose()*m4*J4v+J4a.transpose()*I4*J4a + J5v.transpose()*m5*J5v+J5a.transpose()*I5*J5a + J6v.transpose()*m6*J6v+J6a.transpose()*I6*J6a + Jev.transpose()*me*Jev+Jea.transpose()*Ie*Jea;


    return M;
}

/// do not change the name of the method
inline Eigen::VectorXd getNonlinearities (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {

  /// !!!!!!!!!! NO RAISIM FUNCTIONS HERE !!!!!!!!!!!!!!!!!
  Eigen::VectorXd C, C1, C2, C3, C4, C5, C6, Ce, gravity;
  C1.resize(6);
    C2.resize(6);
    C3.resize(6);
    C4.resize(6);
    C5.resize(6);
    C6.resize(6);
    Ce.resize(6);
    gravity.resize(3);
  C.resize(6);
  gravity << 0, 0, 9.81;
  calculate(gc);
  getJacobian();
  getJacobian_dot(gc, gv);


  C1 << J1v.transpose()*m1*J1v_dot*gv + J1a.transpose()*I1*J1a_dot*gv + J1a.transpose()*CrossProduct((J1a*gv), I1*(J1a*gv)) + J1v.transpose()*gravity*m1;
  C2 << J2v.transpose()*m2*J2v_dot*gv + J2a.transpose()*I2*J2a_dot*gv + J2a.transpose()*CrossProduct((J2a*gv), I2*(J2a*gv)) + J2v.transpose()*gravity*m2;
  C3 << J3v.transpose()*m3*J3v_dot*gv + J3a.transpose()*I3*J3a_dot*gv + J3a.transpose()*CrossProduct((J3a*gv), I3*(J3a*gv)) + J3v.transpose()*gravity*m3;
  C4 << J4v.transpose()*m4*J4v_dot*gv + J4a.transpose()*I4*J4a_dot*gv + J4a.transpose()*CrossProduct((J4a*gv), I4*(J4a*gv)) + J4v.transpose()*gravity*m4;
  C5 << J5v.transpose()*m5*J5v_dot*gv + J5a.transpose()*I5*J5a_dot*gv + J5a.transpose()*CrossProduct((J5a*gv), I5*(J5a*gv)) + J5v.transpose()*gravity*m5;
  C6 << J6v.transpose()*m6*J6v_dot*gv + J6a.transpose()*I6*J6a_dot*gv + J6a.transpose()*CrossProduct((J6a*gv), I6*(J6a*gv)) + J6v.transpose()*gravity*m6;
  Ce << Jev.transpose()*me*Jev_dot*gv + Jea.transpose()*Ie*Jea_dot*gv + Jea.transpose()*CrossProduct((Jea*gv), Ie*(Jea*gv)) + Jev.transpose()*gravity*me;
  C << C1 + C2 + C3 + C4 + C5 + C6 + Ce;

  return C;
}

