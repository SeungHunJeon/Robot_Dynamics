#pragma once

class Contact {
 public:
  Contact() {
    isContact = false;
    lambda_prev.setOnes();
  }
    bool isContact;
    double C_z;
    double C_t;
    Eigen::Vector3d lambda;
    Eigen::Vector3d lambda_prev;
    Eigen::Vector3d lambda_z;
    Eigen::Vector3d lambda_t;
    Eigen::Vector3d v_imp;
    Eigen::Vector3d v_imp_z;
    Eigen::Vector3d v_imp_t;
    Eigen::Vector3d pos;
    Eigen::Matrix3d M_app_inv;
    Eigen::Matrix<double, 3, 6> Jc;
    Eigen::Matrix<double, 6, 1> b;
};

class Ball {
 public:
  Ball() {
    gc.setZero(7);
    gv.setZero(6);
    zero.setZero();
    Identity.setIdentity();
    gv_dot.setZero(6);
    gv_dot << 0, 0, -9.81, 0, 0, 0;
    mu=0.8;
    Jc << Identity, zero;
  }

  double radius;
  double Mass;
  double Inertia;
  double mu;
  Contact ground_contact;
  Contact ball_contact;
  Eigen::Vector3d vel_imp;
  Eigen::VectorXd gc;
  Eigen::VectorXd gv;
  Eigen::VectorXd gv_dot;
  Eigen::Vector3d pos;
  Eigen::Vector3d vel;
  Eigen::Matrix<double, 6, 6> Spatial_mass_matrix;
  Eigen::Matrix<double, 3, 6> Ja;
  Eigen::Matrix<double, 3, 6> Jc;
  Eigen::Matrix<double, 6, 1> b;
  Eigen::Matrix3d zero;
  Eigen::Matrix3d Identity;
  Eigen::Matrix3d Mass_matrix;
  Eigen::Matrix3d Inertia_matrix;
  void cal() {
    Identity_3.setIdentity();
    Zero.setZero();

    Mass = std::pow(radius,3)*M_PI*4/3*rho;
    Inertia = Mass*2/5*std::pow(radius,2);

    Mass_matrix << Identity_3 * Mass;
    Inertia_matrix << Identity_3 * Inertia;

    Spatial_mass_matrix.setZero();
    Spatial_mass_matrix << Mass_matrix, Zero,
    Zero, Inertia_matrix;

    Ja << Zero, Identity_3;
  }

 private:
  double rho=1.0;
  Eigen::Matrix3d Zero;
  Eigen::Matrix3d Identity_3;
};

class SimulationClass {

 public:
  SimulationClass() {
    Identity_3.setIdentity();
    gravity << 0, 0, -9.81;

    ball_1.radius = 0.5;
    ball_1.pos << 0, 0, 0.7;
    ball_1.cal();


    ball_2.radius = 0.7;
    ball_2.pos << 0.1, 0.1, 3;
    ball_2.cal();
  }

  /// hat operator (skew symmetric matrix)
  Eigen::Matrix3d skew(const Eigen::Vector3d vector) {
    Eigen::Matrix3d skew;
    skew << 0, -vector(2), vector(1),
    vector(2), 0, -vector(0),
    -vector(1), vector(0), 0;
    return skew;
  }

  /// world frame into contact frame (get rotation matrix)
  /// Rotation matrix translates the contact normal vector as z-axis
  Eigen::Matrix3d Contact_Rotation(const Eigen::Vector3d axis) {
    Eigen::Matrix3d R, R_z, R_y;
    double theta;
    double phi;
    if(axis(1) >= 0) {
      if(axis.head(2).norm() < 1e-6)
        theta = 0;
      else
        theta = acos(axis(0)/axis.head(2).norm());
    }
    if(axis(1) < 0) {
      if(axis.head(2).norm() < 1e-6)
        theta = 0;
      else
        theta = 2*M_PI - acos(axis(0)/axis.head(2).norm());
    }
    phi = acos(axis(2)/axis.norm());

    R_z << cos(theta), -sin(theta), 0,
    sin(theta), cos(theta), 0,
    0, 0, 1;

    R_y << cos(phi), 0, sin(phi),
    0, 1, 0,
    -sin(phi), 0, cos(phi);

    R << R_z*R_y;

    return R;
  }

  void Contact_integrate() {

    /// no contact on ball_1
    if(Contact_2.isContact == false && Contact_1.isContact == false) {
      ball_1.gv += ball_1.gv_dot * dt;
    }

    /// no contact on ball_2
    if(Contact_3.isContact == false && Contact_1.isContact == false) {
      ball_2.gv += ball_2.gv_dot * dt;
    }

    /// ball_1 - ground contact case
    if(Contact_2.isContact == true && Contact_1.isContact == false) {
      while (true) {
        /// Computation term of v_imp
        ball_1.ground_contact.v_imp << ball_1.gv.head(3) + skew(Contact_2.pos-ball_1.pos)*ball_1.gv.tail(3) + ball_1.ground_contact.M_app_inv*ball_1.ground_contact.lambda -
          ball_1.ground_contact.Jc*ball_1.Spatial_mass_matrix.inverse()*dt*ball_1.ground_contact.b;

        /// set v_imp, v_imp_z, v_imp_t
        ball_1.ground_contact.v_imp_z << 0, 0, (Contact_Rotation(ball_1.pos-Contact_2.pos).transpose()*ball_1.ground_contact.v_imp)(2);
        ball_1.ground_contact.v_imp_t << (Contact_Rotation(ball_1.pos-Contact_2.pos).transpose()*ball_1.ground_contact.v_imp).head(2), 0;
        ball_1.ground_contact.lambda_z << 0, 0, (Contact_Rotation(ball_1.pos-Contact_2.pos).transpose()*ball_1.ground_contact.lambda)(2);
        ball_1.ground_contact.lambda_t << (Contact_Rotation(ball_1.pos-Contact_2.pos).transpose()*ball_1.ground_contact.lambda).head(2), 0;

        /// proj z
        ball_1.ground_contact.lambda_z(2) = std::max(0.0, ball_1.ground_contact.lambda_z(2) - ball_1.ground_contact.v_imp_z(2)*ball_1.ground_contact.C_z);

        /// proj t
        ball_1.ground_contact.lambda_t << ball_1.ground_contact.lambda_t - (ball_1.ground_contact.v_imp_t*ball_1.ground_contact.C_t);
        if(ball_1.ground_contact.lambda_t(0) / ball_1.ground_contact.v_imp_t(0) >= 0)
          ball_1.ground_contact.lambda_t.setZero();
        if(ball_1.ground_contact.lambda_t.norm() >= ball_1.mu * ball_1.ground_contact.lambda_z.norm()) {
          ball_1.ground_contact.lambda_t << ball_1.ground_contact.lambda_t/ball_1.ground_contact.lambda_t.norm()*(ball_1.mu*ball_1.ground_contact.lambda_z.norm());
        }

        /// integrate lambda_t, lambda_z into lambda
        ball_1.ground_contact.lambda << Contact_Rotation(ball_1.pos-Contact_2.pos)*(ball_1.ground_contact.lambda_t + ball_1.ground_contact.lambda_z);

        /// break condition
        if(std::abs((ball_1.ground_contact.lambda-ball_1.ground_contact.lambda_prev).norm()) < 1e-8)
          break;

        /// lambda previous <- lambda
        ball_1.ground_contact.lambda_prev << ball_1.ground_contact.lambda;

      }
      ball_1.gv = ball_1.gv + ball_1.Spatial_mass_matrix.inverse()*dt*(-ball_1.b+ball_1.ground_contact.Jc.transpose()*ball_1.ground_contact.lambda/dt);
    }

    /// ball_2 - ground contact case
    if(Contact_3.isContact == true && Contact_1.isContact == false) {
      while (true) {


        /// Computation term of v_imp
        ball_2.ground_contact.v_imp << ball_2.gv.head(3)+ skew(Contact_3.pos-ball_2.pos)*ball_2.gv.tail(3) + ball_2.ground_contact.M_app_inv*ball_2.ground_contact.lambda -
            ball_2.ground_contact.Jc*ball_2.Spatial_mass_matrix.inverse()*dt*ball_2.ground_contact.b;

        /// set v_imp, v_imp_z, v_imp_t
        ball_2.ground_contact.v_imp_z << 0, 0, (Contact_Rotation(ball_2.pos-Contact_3.pos).transpose()*ball_2.ground_contact.v_imp)(2);
        ball_2.ground_contact.v_imp_t << (Contact_Rotation(ball_2.pos-Contact_3.pos).transpose()*ball_2.ground_contact.v_imp).head(2), 0;
        ball_2.ground_contact.lambda_z << 0, 0, (Contact_Rotation(ball_2.pos-Contact_3.pos).transpose()*ball_2.ground_contact.lambda)(2);
        ball_2.ground_contact.lambda_t << (Contact_Rotation(ball_2.pos-Contact_3.pos).transpose()*ball_2.ground_contact.lambda).head(2), 0;

        /// proj z
        ball_2.ground_contact.lambda_z(2) = std::max(0.0, ball_2.ground_contact.lambda_z(2) - ball_2.ground_contact.v_imp(2)*ball_2.ground_contact.C_z);

        /// proj t
        ball_2.ground_contact.lambda_t << ball_2.ground_contact.lambda_t - (ball_2.ground_contact.v_imp_t*ball_2.ground_contact.C_t);
        if(ball_2.ground_contact.lambda_t(0) / ball_2.ground_contact.v_imp_t(0) >= 0)
          ball_2.ground_contact.lambda_t.setZero();
        if(ball_2.ground_contact.lambda_t.norm() >= ball_2.mu * ball_2.ground_contact.lambda_z.norm()) {
          ball_2.ground_contact.lambda_t << ball_2.ground_contact.lambda_t/ball_2.ground_contact.lambda_t.norm()*(ball_2.mu*ball_2.ground_contact.lambda_z.norm());
        }

        /// integrate lambda_t, lambda_z into lambda
        ball_2.ground_contact.lambda << Contact_Rotation(ball_2.pos-Contact_3.pos)*(ball_2.ground_contact.lambda_t + ball_2.ground_contact.lambda_z);

        /// break condition
        if (std::abs((ball_2.ground_contact.lambda-ball_2.ground_contact.lambda_prev).norm()) < 1e-8)
          break;

        /// lambda previous <- lambda
        ball_2.ground_contact.lambda_prev << ball_2.ground_contact.lambda;

      }
      ball_2.gv = ball_2.gv + ball_2.Spatial_mass_matrix.inverse()*dt*(-ball_2.b+ball_2.ground_contact.Jc.transpose()*ball_2.ground_contact.lambda/dt);


    }

    /// One of the specified case that ball - ball contact && ball 1 - ground contact Simultaeneously
    if(Contact_1.isContact == true && Contact_2.isContact == true) {

      while (true) {
        /// Computation term of v_imp
        ball_1.ground_contact.v_imp << ball_1.gv.head(3) + skew(Contact_2.pos-ball_1.pos)*ball_1.gv.tail(3)
        + ball_1.ground_contact.M_app_inv*ball_1.ground_contact.lambda
        + ball_1.ground_contact.Jc*ball_1.Spatial_mass_matrix.inverse()*ball_1.ball_contact.Jc.transpose()*ball_1.ball_contact.lambda
        - dt*ball_1.ground_contact.Jc*ball_1.Spatial_mass_matrix.inverse()*ball_1.ground_contact.b;

        ball_1.ball_contact.v_imp << (ball_1.gv.head(3) + skew(Contact_1.pos-ball_1.pos)*ball_1.gv.tail(3))-(ball_2.gv.head(3) + skew(Contact_1.pos-ball_2.pos)*ball_2.gv.tail(3))
        + (ball_1.ball_contact.M_app_inv + ball_2.ball_contact.M_app_inv)*ball_1.ball_contact.lambda
        + ball_1.ball_contact.Jc*ball_1.Spatial_mass_matrix.inverse()*ball_1.ground_contact.Jc.transpose()*ball_1.ground_contact.lambda
        + dt*(-ball_1.ball_contact.Jc*ball_1.Spatial_mass_matrix.inverse()*ball_1.ball_contact.b + ball_2.ball_contact.Jc*ball_2.Spatial_mass_matrix.inverse()*ball_2.ball_contact.b);

        /// set v_imp, v_imp_t, v_imp_z
        ball_1.ball_contact.v_imp_z << 0, 0, (Contact_Rotation(ball_1.pos-Contact_1.pos).transpose()*ball_1.ball_contact.v_imp)(2);
        ball_1.ball_contact.v_imp_t << (Contact_Rotation(ball_1.pos-Contact_1.pos).transpose()*ball_1.ball_contact.v_imp).head(2), 0;
        ball_1.ball_contact.lambda_z << 0, 0, (Contact_Rotation(ball_1.pos-Contact_1.pos).transpose()*ball_1.ball_contact.lambda)(2);
        ball_1.ball_contact.lambda_t << (Contact_Rotation(ball_1.pos-Contact_1.pos).transpose()*ball_1.ball_contact.lambda).head(2), 0;

        ball_1.ground_contact.v_imp_z << 0, 0, (Contact_Rotation(ball_1.pos-Contact_2.pos).transpose()*ball_1.ground_contact.v_imp)(2);
        ball_1.ground_contact.v_imp_t << (Contact_Rotation(ball_1.pos-Contact_2.pos).transpose()*ball_1.ground_contact.v_imp).head(2), 0;
        ball_1.ground_contact.lambda_z << 0, 0, (Contact_Rotation(ball_1.pos-Contact_2.pos).transpose()*ball_1.ground_contact.lambda)(2);
        ball_1.ground_contact.lambda_t << (Contact_Rotation(ball_1.pos-Contact_2.pos).transpose()*ball_1.ground_contact.lambda).head(2), 0;


        /// proj z
        ball_1.ball_contact.lambda_z(2) = std::max(0.0, ball_1.ball_contact.lambda_z(2) - ball_1.ball_contact.v_imp_z(2)*ball_1.ball_contact.C_z);
        ball_1.ground_contact.lambda_z(2) = std::max(0.0, ball_1.ground_contact.lambda_z(2) - ball_1.ground_contact.v_imp_z(2)*ball_1.ground_contact.C_z);

        /// proj t
        ball_1.ground_contact.lambda_t << ball_1.ground_contact.lambda_t - (ball_1.ground_contact.v_imp_t*ball_1.ground_contact.C_t);
        if(ball_1.ground_contact.lambda_t(0) / ball_1.ground_contact.v_imp_t(0) >= 0)
          ball_1.ground_contact.lambda_t.setZero();
        if(ball_1.ground_contact.lambda_t.norm() >= ball_1.mu * ball_1.ground_contact.lambda_z.norm()) {
          ball_1.ground_contact.lambda_t << ball_1.ground_contact.lambda_t/ball_1.ground_contact.lambda_t.norm()*(ball_1.mu*ball_1.ground_contact.lambda_z.norm());
        }

        ball_1.ball_contact.lambda_t << ball_1.ball_contact.lambda_t - (ball_1.ball_contact.v_imp_t*ball_1.ball_contact.C_t);
        if(ball_1.ball_contact.lambda_t(0) / ball_1.ball_contact.v_imp_t(0) >= 0)
          ball_1.ball_contact.lambda_t.setZero();

        if(ball_1.ball_contact.lambda_t.norm() >= ball_1.mu * ball_1.ball_contact.lambda_z.norm()) {
            ball_1.ball_contact.lambda_t << ball_1.ball_contact.lambda_t/(ball_1.ball_contact.lambda_t.norm())*(ball_1.mu * ball_1.ball_contact.lambda_z.norm());
        }

        /// integrate lambda with lambda_t & lambda_z
        ball_1.ball_contact.lambda << (Contact_Rotation(ball_1.pos-Contact_1.pos)*(ball_1.ball_contact.lambda_t + ball_1.ball_contact.lambda_z));
        ball_1.ground_contact.lambda << Contact_Rotation(ball_1.pos-Contact_2.pos)*(ball_1.ground_contact.lambda_t + ball_1.ground_contact.lambda_z);
        ball_2.ball_contact.lambda << -ball_1.ball_contact.lambda;

        /// break condition
        if((std::abs((ball_1.ball_contact.lambda-ball_1.ball_contact.lambda_prev).norm()) + std::abs((ball_1.ground_contact.lambda-ball_1.ground_contact.lambda_prev).norm())) < 1e-8)
          break;

        /// lambda previous <- lambda
        ball_1.ball_contact.lambda_prev << ball_1.ball_contact.lambda;
        ball_1.ground_contact.lambda_prev << ball_1.ground_contact.lambda;
      }

      ball_1.gv = ball_1.gv + ball_1.Spatial_mass_matrix.inverse()*dt*(-ball_1.b+ball_1.ball_contact.Jc.transpose()*ball_1.ball_contact.lambda/dt+ball_1.ground_contact.Jc.transpose()*ball_1.ground_contact.lambda/dt);
      ball_2.gv = ball_2.gv + ball_2.Spatial_mass_matrix.inverse()*dt*(-ball_2.b+ball_2.ball_contact.Jc.transpose()*ball_2.ball_contact.lambda/dt);
    }





  };

  void isContact() {
    /// contact detection, contact_1 : ball to ball, contact_2 : ball_1 to ground, contact_3 : ball_2 to ground

    if ((ball_1.pos - ball_2.pos).norm() < (ball_1.radius + ball_2.radius)) {
      Contact_1.isContact = true;

      /// contact pos
      Contact_1.pos << ball_1.pos + ball_1.radius * (ball_2.pos - ball_1.pos) / ((ball_2.pos - ball_1.pos).norm());

      /// compute ball_1
      ball_1.ball_contact.Jc << Identity_3, skew(Contact_1.pos-ball_1.pos);
      ball_1.ball_contact.M_app_inv << ball_1.ball_contact.Jc*ball_1.Spatial_mass_matrix.inverse()*ball_1.ball_contact.Jc.transpose();
      ball_1.ball_contact.C_z = constant/ball_1.ball_contact.M_app_inv(8);
      ball_1.ball_contact.C_t = constant/std::max(ball_1.ball_contact.M_app_inv(0),ball_1.ball_contact.M_app_inv(4));
      ball_1.ball_contact.b << ball_1.Ja.transpose()*skew(ball_1.gv.tail(3))*(ball_1.Inertia_matrix*ball_1.gv.tail(3)) - ball_1.Jc.transpose()*gravity*ball_1.Mass;
      ball_1.b << ball_1.Ja.transpose()*skew(ball_1.gv.tail(3))*(ball_1.Inertia_matrix*ball_1.gv.tail(3)) - ball_1.Jc.transpose()*gravity*ball_1.Mass;
      /// compute ball_2
      ball_2.ball_contact.Jc << Identity_3, skew(Contact_1.pos-ball_2.pos);
      ball_2.ball_contact.M_app_inv << ball_2.ball_contact.Jc*ball_2.Spatial_mass_matrix.inverse()*ball_2.ball_contact.Jc.transpose();
      ball_2.ball_contact.C_z = constant/ball_2.ball_contact.M_app_inv(8);
      ball_2.ball_contact.C_t = constant/std::max(ball_2.ball_contact.M_app_inv(0),ball_2.ball_contact.M_app_inv(4));
      ball_2.ball_contact.b << ball_2.Ja.transpose()*skew(ball_2.gv.tail(3))*(ball_2.Inertia_matrix*ball_2.gv.tail(3)) - ball_2.Jc.transpose()*gravity*ball_2.Mass;
      ball_2.b << ball_2.Ja.transpose()*skew(ball_2.gv.tail(3))*(ball_2.Inertia_matrix*ball_2.gv.tail(3)) - ball_2.Jc.transpose()*gravity*ball_2.Mass;
    } else
      Contact_1.isContact = false;
    if (ball_1.pos(2) < ball_1.radius) {
      Contact_2.isContact = true;
      /// compute Contact 2 pos
      Contact_2.pos << ball_1.pos(0), ball_1.pos(1), 0;
      /// compute ball_1
      ball_1.ground_contact.Jc  << Identity_3, skew(Contact_2.pos-ball_1.pos);
      ball_1.ground_contact.M_app_inv << ball_1.ground_contact.Jc*ball_1.Spatial_mass_matrix.inverse()*ball_1.ground_contact.Jc.transpose();
      ball_1.ground_contact.C_z = constant/ball_1.ground_contact.M_app_inv(8);
      ball_1.ground_contact.C_t = constant/std::max(ball_1.ground_contact.M_app_inv(0),ball_1.ground_contact.M_app_inv(4));
      ball_1.ground_contact.b << ball_1.Ja.transpose()*skew(ball_1.gv.tail(3))*(ball_1.Inertia_matrix*ball_1.gv.tail(3)) - ball_1.Jc.transpose()*gravity*ball_1.Mass;
      ball_1.b << ball_1.Ja.transpose()*skew(ball_1.gv.tail(3))*(ball_1.Inertia_matrix*ball_1.gv.tail(3)) - ball_1.Jc.transpose()*gravity*ball_1.Mass;

    } else
      Contact_2.isContact = false;

    if (ball_2.pos(2) < ball_2.radius) {
      Contact_3.isContact = true;
      /// compute Contact 3 pos
      Contact_3.pos << ball_2.pos(0), ball_2.pos(1), 0;
      /// compute ball_2
      ball_2.ground_contact.Jc << Identity_3, skew(Contact_3.pos-ball_2.pos);
      ball_2.ground_contact.M_app_inv << ball_2.ground_contact.Jc*ball_2.Spatial_mass_matrix.inverse()*ball_2.ground_contact.Jc.transpose();
      ball_2.ground_contact.C_z = constant/ball_2.ground_contact.M_app_inv(8);
      ball_2.ground_contact.C_t = constant/std::max(ball_2.ground_contact.M_app_inv(0),ball_2.ground_contact.M_app_inv(4));
      ball_2.ground_contact.b << ball_2.Ja.transpose()*skew(ball_2.gv.tail(3))*(ball_2.Inertia_matrix*ball_2.gv.tail(3)) - ball_2.Jc.transpose()*gravity*ball_2.Mass;
      ball_2.b << ball_2.Ja.transpose()*skew(ball_2.gv.tail(3))*(ball_2.Inertia_matrix*ball_2.gv.tail(3)) - ball_2.Jc.transpose()*gravity*ball_2.Mass;

    } else
      Contact_3.isContact = false;
  }

  void integrate() {
    isContact();
    Contact_integrate();
    ball_1.pos << ball_1.pos + ball_1.gv.segment(0, 3) * dt;
    ball_2.pos << ball_2.pos + ball_2.gv.segment(0,3)*dt;
  }

  void setPosition(raisim::Visuals * sphere1, raisim::Visuals * sphere2) {
    sphere1->setPosition(ball_1.pos);
    sphere2->setPosition(ball_2.pos);
  }
 private:
  double dt = 0.001;
  double constant = 0.01;
  Eigen::Vector3d gravity;
  Ball ball_1, ball_2;
  Contact Contact_1, Contact_2, Contact_3;
  Eigen::Matrix3d Identity_3;
  Eigen::Matrix3d Zero;

  /// add state variables here
};