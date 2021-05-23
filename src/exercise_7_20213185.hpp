#pragma once


class SimulationClass {

 public:
  SimulationClass(){
    ball.resize(2);
    ball[0].radius=0.5;
    ball[1].radius=0.7;
    ball[0].pos_={0,0,0.5};
    ball[1].pos_={0.1,0.1,3};
    ball[0].pos__={0,0,0};
    ball[1].pos__={0,0,0};
    ball[0].vel_={0,0,0};
    ball[1].vel_={0,0,0};
    ball[0].vel__={0,0,0};
    ball[1].vel__={0,0,0};
    ball[0].acc=-g;
    ball[1].acc=-g;
    I_mat_6.setIdentity(6,6);
    Ja.setZero(3,6);
    Ja.block(0,3,3,3)=I_mat_3;
    for (int i=0;i<2;i++){
      ball[i].mass=4*M_PI/3.*pow(ball[i].radius,3);
      ball[i].w= {0,0,0};
      ball[i].Jg.setZero(3,6);
      ball[i].Jb.setZero(3,6);
      ball[i].R.setIdentity(3,3);
      ball[i].M=ball[i].mass*I_mat_6;
      ball[i].M.block(3,3,3,3)=2/5.*ball[i].mass*pow(ball[i].radius,2)*I_mat_3;
      ball[i].lambda_ball.setZero();
      ball[i].lambda_ground.setZero();
      ball[i].contact_point_b.setZero();
    }
  } // 생성자
  static Eigen::Matrix3d hat(Eigen::Vector3d v){
    Eigen::Matrix3d vv;
    vv << 0, -v(2), v(1),
            v(2), 0, -v(0),
            -v(1), v(0),0;
    return vv;
  }
  void integrate() {
    freefall(0);//dt동안 자유낙하
    freefall(1);//dt동안 자유낙하

//    for (int i=0;i<2;i++){
//      std::cout<<"Contatct mode of 0: "<<Contact_mode(0)<<std::endl;
//      switch (Contact_mode(i)){
//        case 0:
//          break;
//        case 1:
//          Collideground(i);
//          break;
//        case 2:
//          if (i==0) Collideball(0,1);
//          else Collideball(0,1);
//          break;
//        default:
//          break;
//      }
//      ball[i].pos_=ball[i].pos__;
//      ball[i].vel_=ball[i].vel__;
//    }
    std::cout<<"test_case_num:"<<test_case_num()<<std::endl;
    switch(test_case_num()){
      case 1:
        Collideground(0);
        std::cout<<"case 1"<<std::endl;
        break;
      case 2:
        CollideCase2();
        break;
      case 3:
        CollideCase3();
        break;
      case 4:
        Collideground(0);
        Collideground(1);
        break;
      case 5:
        break;
      case 6:
        Collideball(0,1);
        break;
      case 7:
        Collideground(1);
        break;
    }

    for (int i=0;i<2;i++) {
      ball[i].pos_ = ball[i].pos__;
      ball[i].vel_ = ball[i].vel__;
    }
  }
  int test_case_num(){ //오직 지금 config에만 적용
    if(Contact_mode(0)==1 and Contact_mode(1)==0) return 1;
    if(Contact_mode(0)==3 and Contact_mode(1)==2) return 2;
    if(Contact_mode(0)==3 and Contact_mode(1)==3) return 3;
    if(Contact_mode(0)==1 and Contact_mode(1)==1) return 4;
    if(Contact_mode(0)==0 and Contact_mode(1)==0) return 5;
    if(Contact_mode(0)==2) return 6;
    if(Contact_mode(0)==0 and Contact_mode(1)==1) return 7;
    else {
      std::cout<<"test_case_num_error"<<std::endl;
      std::cout<<"Contact_mode of ball0: "<<Contact_mode(0)<<std::endl;
      std::cout<<"Contact_mode of ball1: "<<Contact_mode(1)<<std::endl;
      return 100;
    }
  }
  void Collideground(int i){
    Eigen::Vector3d Vc_before,Vc_after, rc, VB;
    Eigen::VectorXd bterm(6);
    double lx=0, lz=0,lx_prev=0,lz_prev=0,prox_x=0,prox_z=0, alpha=0.1, beta=0.1;
    double error_th=1e-8;
    rc<< 0,0,-ball[i].radius;

    Vc_before=ball[i].vel_+hat(ball[i].w)*rc;
    VB<<Vc_before(0),Vc_before(1),0;
    if (VB.norm()>1e-6) {
      ball[i].R.col(0) = VB / VB.norm(); // 충돌직전 공의 지면에 대한 tangential상대속도의 방향을 x축으로 잡는다.
      ball[i].R.col(2) << 0, 0, 1; //ground 에서 contact의 z축은 world의 z축
      ball[i].R.col(1) = hat(ball[i].R.col(2)) * ball[i].R.col(0); // y축은 cross(z,x)
      ball[i].R.col(1) = (ball[i].R.col(1)) / (ball[i].R.col(1)).norm();
    }
    ball[i].Jg.block(0,0,3,3)=I_mat_3;
    ball[i].Jg.block(0,3,3,3)=-hat(rc);
    bterm<<ball[i].mass*g,0,0,0;
    bterm+=Ja.transpose()*(hat(ball[i].w)*(ball[i].M.block(3,3,3,3)*ball[i].w));
    while (true){
      ball[i].lambda_ground << lx,0,lz; //contact frame기준
      Vc_after=Vc_before+(ball[i].Jg*(ball[i].M.inverse())*ball[i].Jg.transpose())*(ball[i].R*ball[i].lambda_ground)-dt*(ball[i].Jg*ball[i].M.inverse()*bterm);
      prox_z=lz-alpha*(ball[i].R.transpose()*Vc_after)(2)*(1/ball[i].M(2,2)); //contact frame기준
      if (prox_z<0) prox_z=0;
      lz=prox_z;
      ball[i].lambda_ground << lx,0,lz;
      Vc_after=Vc_before+(ball[i].Jg*(ball[i].M.inverse())*ball[i].Jg.transpose())*(ball[i].R*ball[i].lambda_ground)-dt*(ball[i].Jg*ball[i].M.inverse()*bterm);
      prox_x=lx-beta*(ball[i].R.transpose()*Vc_after)(0)/(fmax(1/ball[i].M(1,1),1/ball[i].M(2,2)));
      if (prox_x>0) prox_x=0;
      if (prox_x<-mu*lz) prox_x=-mu*lz;
      lx=prox_x;
      if(pow(lx-lx_prev,2)<error_th and pow(lz-lz_prev,2)<error_th) { //충분히 수렴한 경우
        ball[i].lambda_ground << lx,0,lz;
        ball[i].vel__=ball[i].vel_+(ball[i].M.inverse()*((-bterm)*dt+ball[i].Jg.transpose()*ball[i].R*ball[i].lambda_ground)).block(0,0,3,1);
        ball[i].w+=(ball[i].M.inverse()*((-bterm)*dt+ball[i].Jg.transpose()*ball[i].R*ball[i].lambda_ground)).block(3,0,3,1);
        ball[i].pos__=ball[i].pos_+dt*ball[i].vel__;

        break;
      }
      lx_prev=lx;
      lz_prev=lz;
    }

  }
  double innerproduct(Eigen::Vector3d v1, Eigen::Vector3d v2){
    double sum=0;
    for (int i=0;i<3;i++){
      sum+=v1(i)*v2(i);
    }
    return sum;
  }
  void Collideball(int i, int j) { // ball1이 i ,ball2가 j인 경우만 생각하자.
    Eigen::Vector3d Vc_before, Vc_after, rc_i,rc_j,VB,rc_i_normalized;
    Eigen::VectorXd bterm_i(6), bterm_j(6), bterm_ij(6);
    Eigen::Matrix3d Mapp;
    double lx = 0, lz = 0, lx_prev = 0, lz_prev = 0, prox_x = 0, prox_z = 0, alpha = 0.5, beta = 0.5;
    double error_th = 1e-10;
    rc_i << (ball[j].pos_-ball[i].pos_)/(ball[i].radius+ball[j].radius)*ball[i].radius;
    rc_j << (ball[i].pos_-ball[j].pos_)/(ball[j].radius+ball[i].radius)*ball[j].radius;
    rc_i_normalized = rc_i/rc_i.norm();
    Vc_before = (ball[i].vel_ + hat(ball[i].w) * rc_i)-(ball[j].vel_ + hat(ball[j].w) * rc_j); //Vc_before은 world frame기준

    VB = Vc_before- innerproduct(Vc_before,rc_i_normalized)*rc_i_normalized;
    if (VB.norm() > 1e-8) {
      ball[i].R.col(0) = VB / VB.norm(); // 충돌직전 공의 지면에 대한 tangential상대속도의 방향을 x축으로 잡는다.
      ball[i].R.col(2) =-rc_i_normalized; //ground 에서 contact의 z축은 world의 z축
      ball[i].R.col(1) = hat(ball[i].R.col(2)) * ball[i].R.col(0); // y축은 cross(z,x)
      ball[i].R.col(1) = (ball[i].R.col(1)) / (ball[i].R.col(1)).norm();
      ball[j].R=-ball[i].R;
    }

    ball[i].Jg.block(0, 0, 3, 3) = I_mat_3;
    ball[i].Jg.block(0, 3, 3, 3) = -hat(rc_i);
    ball[j].Jg.block(0, 0, 3, 3) = I_mat_3;
    ball[j].Jg.block(0, 3, 3, 3) = -hat(rc_j);

    bterm_i << ball[i].mass * g,0,0,0;
    bterm_i +=Ja.transpose() * (hat(ball[i].w) * (ball[i].M.block(3, 3, 3, 3) * ball[i].w));
    bterm_j << ball[j].mass * g,0,0,0;
    bterm_j +=Ja.transpose() * (hat(ball[j].w) * (ball[j].M.block(3, 3, 3, 3) * ball[j].w));

    Mapp=((ball[i].Jg * (ball[i].M.inverse()) * ball[i].Jg.transpose())+((ball[j].Jg * (ball[j].M.inverse()) * ball[j].Jg.transpose()))).inverse();

    while (true) {
      //std::cout<<"in while loop"<<std::endl;
      ball[i].lambda_ball << lx, 0, lz; //contact frame기준

      Vc_after = Vc_before + Mapp.inverse()*ball[i].R*ball[i].lambda_ball
                -dt * (ball[i].Jg * ball[i].M.inverse() * bterm_i)
                +dt * (ball[j].Jg * ball[j].M.inverse() * bterm_j);
      prox_z = lz - alpha * (ball[i].R.transpose() * Vc_after)(2) / (1 / Mapp(2, 2)); //contact frame기준
      if (prox_z < 0) prox_z = 0;
      lz = prox_z;
      ball[i].lambda_ball << lx, 0, lz; //lambda ball은 contact frame기준
      Vc_after = Vc_before + Mapp.inverse()*ball[i].R*ball[i].lambda_ball
                 -dt * (ball[i].Jg * ball[i].M.inverse() * bterm_i)
                 +dt * (ball[j].Jg * ball[j].M.inverse() * bterm_j);
      prox_x = lx - beta * (ball[i].R.transpose() * Vc_after)(0) / (fmax(1 / Mapp(0, 0), 1 / Mapp(1, 1)));

      if (prox_x > 0) prox_x = 0;
      if (prox_x < -mu * lz) prox_x = -mu * lz;
      lx = prox_x;
      if (pow(lx - lx_prev, 2) < error_th and pow(lz - lz_prev, 2) < error_th) { //충분히 수렴한 경우

        ball[i].lambda_ball << lx, 0, lz;
        ball[j].lambda_ball = ball[i].lambda_ball;
        ball[i].vel__ = ball[i].vel_ + (ball[i].M.inverse() * ((-bterm_i) * dt + ball[i].Jg.transpose() * ball[i].R *
                                                                               ball[i].lambda_ball)).block(0, 0, 3,
                                                                                                             1);
        ball[i].w += (ball[i].M.inverse() *
                      ((-bterm_i) * dt + ball[i].Jg.transpose() * ball[i].R * ball[i].lambda_ball)).block(3, 0, 3, 1);
        ball[i].pos__ = ball[i].pos_ + dt * ball[i].vel__;


        ball[j].vel__ = ball[j].vel_ + (ball[j].M.inverse() * ((-bterm_j) * dt + ball[j].Jg.transpose() * ball[j].R *
                                                                                 ball[j].lambda_ball)).block(0, 0, 3,
                                                                                                               1);
        ball[j].w += (ball[j].M.inverse() *
                      ((-bterm_j) * dt + ball[j].Jg.transpose() * ball[j].R * ball[j].lambda_ball)).block(3, 0, 3, 1);
        ball[j].pos__ = ball[j].pos_ + dt * ball[j].vel__;

        break;
      }
      lx_prev = lx;
      lz_prev = lz;
    }
  }
  void CollideCase2(){
    Eigen::Vector3d Vc_before_1,Vc_before_2, Vc_after_1,Vc_after_2, rc_01,rc_02, rc_12, VB,VB_, rc_02_normalized;
    Eigen::VectorXd bterm_0(6), bterm_1(6);
    Eigen::Matrix3d Mapp1, Mapp02, Mapp12, Mapp2_inv;
    Eigen::Vector3d lambda1, lambda2;
    Eigen::Matrix3d R01,R02,R12;
    Eigen::VectorXd lambda_now(4),lambda_prev(4);
    R01.setIdentity();
    R02.setIdentity();
    R12.setIdentity();

    double lx1=0,lx2=0,lz1=0,lz2=0,lx_prev1=0,lx_prev2=0,lz_prev1=0,lz_prev2=0;
    double prox_x1 = 0,prox_x2 = 0,prox_z1 = 0, prox_z2 = 0;
    double alpha = 0.5, beta = 0.5;
    double error_th = 1e-10;
    rc_01={0,0,-ball[0].radius};
    rc_02 << (ball[1].pos_ - ball[0].pos_) / (ball[0].radius + ball[1].radius) * ball[0].radius;
    rc_12 << (ball[0].pos_ - ball[1].pos_) / (ball[0].radius + ball[1].radius) * ball[1].radius;
//    std::cout<<"rc_01: "<<rc_01<<std::endl;
//    std::cout<<"rc_02: "<<rc_02<<std::endl;
//    std::cout<<"rc_12: "<<rc_12<<std::endl;

    rc_02_normalized = rc_02 / rc_02.norm();
    Vc_before_1 =(ball[0].vel_ + hat(ball[0].w) * rc_01); //Vc_before은 world frame기준
    Vc_before_2=(ball[0].vel_ + hat(ball[0].w) * rc_02)-(ball[1].vel_ + hat(ball[1].w) * rc_12);


    VB = Vc_before_2 - innerproduct(Vc_before_2, rc_02_normalized) * rc_02_normalized;
    std::cout<<"VB: "<<VB<<std::endl;

    if (VB.norm() > 1e-8) {
      R02.col(0) = VB / VB.norm(); // 충돌직전 공의 지면에 대한 tangential상대속도의 방향을 x축으로 잡는다.
      R02.col(2) = -rc_02_normalized; //ground 에서 contact의 z축은 world의 z축
      R02.col(1) = hat(R02.col(2)) * R02.col(0); // y축은 cross(z,x)
      R02.col(1) = (R02.col(1)) / (R02.col(1)).norm();
      R12=-R02;
    }

    VB_<< Vc_before_1(0),Vc_before_1(1),0; //TODO
    if (VB_.norm()>1e-8) {
      R01.col(0) = Vc_before_1 / Vc_before_1.norm(); // 충돌직전 공의 지면에 대한 tangential상대속도의 방향을 x축으로 잡는다.
      R01.col(2) << 0, 0, 1; //ground 에서 contact의 z축은 world의 z축
      R01.col(1) = hat(R01.col(2)) * R01.col(0); // y축은 cross(z,x)
      R01.col(1) = (R01.col(1)) / (R01.col(1)).norm();
    }

    std::cout<<R01<<std::endl;
    std::cout<<R02<<std::endl;
    std::cout<<Vc_before_2<<std::endl;

    ball[0].Jg.block(0, 0, 3, 3) = I_mat_3;
    ball[0].Jg.block(0, 3, 3, 3) = -hat(rc_01);
    ball[0].Jb.block(0, 0, 3, 3) = I_mat_3;
    ball[0].Jb.block(0, 3, 3, 3) = -hat(rc_02);
    ball[1].Jb.block(0, 0, 3, 3) = I_mat_3;
    ball[1].Jb.block(0, 3, 3, 3) = -hat(rc_12);

    bterm_0 << ball[0].mass * g, 0, 0, 0;
    bterm_0 += Ja.transpose() * (hat(ball[0].w) * (ball[0].M.block(3, 3, 3, 3) * ball[0].w));
    bterm_1 << ball[1].mass * g, 0, 0, 0;
    bterm_1 += Ja.transpose() * (hat(ball[1].w) * (ball[1].M.block(3, 3, 3, 3) * ball[1].w));

    Mapp1= (ball[0].Jg * (ball[0].M.inverse()) * ball[0].Jg.transpose()).inverse();
    Mapp02 = (ball[0].Jb * (ball[0].M.inverse()) * ball[0].Jb.transpose()).inverse();
    Mapp12 = (ball[1].Jb * (ball[1].M.inverse()) * ball[1].Jb.transpose()).inverse();
    Mapp2_inv=Mapp02.inverse()+Mapp12.inverse();
//    std::cout<<"Vc_before_1: "<<Vc_before_1<<std::endl;
//    std::cout<<"Vc_before_2: "<<Vc_before_2<<std::endl;
    while (true) { ///////////////////////////////////////////////////////////////////////////////////////////////////////////

      lambda1 << lx1, 0, lz1; //contact frame기준
      lambda2 << lx2, 0, lz2; //contact frame기준


      Vc_after_1= Vc_before_1 + (ball[0].Jg*(ball[0].M.inverse())*ball[0].Jg.transpose())*(R01*lambda1)
                  +(ball[0].Jg*(ball[0].M.inverse())*ball[0].Jb.transpose())*(R02*lambda2)
                  -dt*(ball[0].Jg*ball[0].M.inverse()*bterm_0);
      prox_z1 = lz1 - alpha * (R01.transpose() * Vc_after_1)(2) / (1 / Mapp1(2, 2)); //contact frame기준
      prox_x1 = lx1 - beta * (R01.transpose() * Vc_after_1)(0) / fmax(1 / Mapp1(0,0),1 / Mapp1(1, 1)); //contact frame기준
      lz1= cons_prox_z(prox_z1);
      lx1= cons_prox_x(prox_x1,lz1);
//      std::cout<<"-dt*(ball[0].Jg*ball[0].M.inverse()*bterm_0): "<<-dt*(ball[0].Jg*ball[0].M.inverse()*bterm_0)<<std::endl;
//      std::cout<<"proxx1: "<<prox_x1<<std::endl;
//      std::cout<<"Vc_after_1: "<<Vc_after_1<<std::endl;

      lambda1 << lx1, 0, lz1; //contact frame기준


      Vc_after_2= Vc_before_2 +
              (ball[0].Jb*(ball[0].M.inverse())*ball[0].Jg.transpose())*(R01*lambda1)
                  +((ball[0].Jb * (ball[0].M.inverse()) * ball[0].Jb.transpose())+(ball[1].Jb * (ball[1].M.inverse()) * ball[1].Jb.transpose()))*(R02*lambda2)
                  -dt * (ball[0].Jb * ball[0].M.inverse() * bterm_0)
                  +dt * (ball[1].Jb * ball[1].M.inverse() * bterm_1);



      prox_x2 = lx2 - beta * (R02.transpose() * Vc_after_2)(0) / fmax(Mapp2_inv(0,0),Mapp2_inv(1, 1)); //contact frame기준
      prox_z2 = lz2 - alpha * (R02.transpose() * Vc_after_2)(2) / (Mapp2_inv(2, 2)); //contact frame기준
      lz2= cons_prox_z(prox_z2);
      lx2= cons_prox_x(prox_x2,lz2);

      lambda2 << lx2, 0, lz2; //contact frame기준


      lambda_prev<<lx_prev1,lx_prev2,lz_prev1,lz_prev2;
      lambda_now<<lx1,lx2,lz1,lz2;


      if ((lambda_prev-lambda_now).norm() < error_th) { //충분히 수렴한 경우
//        std::cout<<"Vc_after1: "<<Vc_after_1<<std::endl;
//        std::cout<<"Vc_after2: "<<Vc_after_2<<std::endl;
//        std::cout<<"in :"<<lambda_now<<std::endl;
//        std::cout<<"lambda_now:"<<lambda_now<<std::endl;
        lambda1 << lx1, 0, lz1; //contact frame기준
        lambda2 << lx2, 0, lz2; //contact frame기준

//        std::cout<<"Vc_after1: "<<Vc_after_1<<std::endl;
//        std::cout<<"Vc_after2: "<<Vc_after_2<<std::endl;

        ball[0].vel__ = ball[0].vel_ + (ball[0].M.inverse() * ((-bterm_0) * dt
                                                               + ball[0].Jg.transpose() * R01 *lambda1+ball[0].Jb.transpose() * R02 *lambda2)).block(0, 0, 3,1);
        ball[0].w += (ball[0].M.inverse() * ((-bterm_0) * dt
                                             + ball[0].Jg.transpose() * R01 *lambda1+ball[0].Jb.transpose() * R02 *lambda2)).block(3, 0, 3,1);
        ball[0].pos__ = ball[0].pos_ + dt * ball[0].vel__;

        ball[1].vel__ = ball[1].vel_ + (ball[1].M.inverse() * ((-bterm_1) * dt + ball[1].Jb.transpose() * R12 *
                                                                                                 lambda2)).block(0, 0, 3,
                                                                                                                             1);
        ball[1].w += (ball[1].M.inverse() * ((-bterm_1) * dt
                                             +ball[1].Jb.transpose() * R12 *lambda2)).block(3, 0, 3,1);
        ball[1].pos__ = ball[1].pos_ + dt * ball[1].vel__;

//        std::cout<<"ball[0].vel__: "<<ball[0].vel__<<std::endl;
//        std::cout<<"ball[1].vel__: "<<ball[1].vel__<<std::endl;

        break;
      }
      lx_prev1 = lx1;
      lz_prev1 = lz1;
      lx_prev2 = lx2;
      lz_prev2 = lz2;
    }
  }
  double cons_prox_z(double z){
    if(z<0){
      z=0;
    }
    return z;
  }
  double cons_prox_x(double x, double z){
    if (x > 0){
      x=0;
    }
    if(x<-mu*z){
      x=-mu*z;
    }
    return x;
  }
  void CollideCase3() {
    int i=0, j=1;
    Eigen::Vector3d Vc_before_1,Vc_before_2,Vc_before_3, Vc_after_1,Vc_after_2,Vc_after_3, rc_01,rc_02, rc_12, rc_13, VB, rc_02_normalized;
    Eigen::VectorXd bterm_0(6), bterm_1(6);
    Eigen::Matrix3d Mapp;
    Eigen::Vector3d lambda1, lambda2, lambda3;
    Eigen::Matrix3d R01,R02,R12,R13;
    Eigen::VectorXd lambda_now(6),lambda_prev(6);
    R01.setIdentity(),R02.setIdentity(),R12.setIdentity(),R13.setIdentity();

    double lx1=0,lx2=0,lx3=0,lz1=0,lz2=0,lz3=0,lx_prev1=0,lx_prev2=0,lx_prev3=0,lz_prev1=0,lz_prev2=0,lz_prev3=0;
    double prox_x1 = 0,prox_x2 = 0,prox_x3 = 0,prox_z1 = 0, prox_z2 = 0,prox_z3 = 0;
    double alpha = 0.1, beta = 0.1;
    double error_th = 1e-8;
    rc_01={0,0,-ball[0].radius};
    rc_13={0,0,-ball[1].radius};
    rc_02 << (ball[1].pos_ - ball[0].pos_) / (ball[0].radius + ball[1].radius) * ball[0].radius;
    rc_12 << (ball[0].pos_ - ball[1].pos_) / (ball[0].radius + ball[1].radius) * ball[1].radius;
    rc_02_normalized = rc_02 / rc_02.norm();
    Vc_before_1 =(ball[0].vel_ + hat(ball[0].w) * rc_01); //Vc_before은 world frame기준
    Vc_before_2=(ball[0].vel_ + hat(ball[0].w) * rc_01)-(ball[1].vel_ + hat(ball[1].w) * rc_12);
    Vc_before_3=(ball[1].vel_ + hat(ball[1].w) * rc_13);

    VB = Vc_before_2 - innerproduct(Vc_before_2, rc_02_normalized) * rc_02_normalized;

    if (VB.norm() > 1e-8) {
      R02.col(0) = VB / VB.norm(); // 충돌직전 공의 지면에 대한 tangential상대속도의 방향을 x축으로 잡는다.
      R02.col(2) = -rc_02_normalized; //ground 에서 contact의 z축은 world의 z축
      R02.col(1) = hat(R02.col(2)) * R02.col(0); // y축은 cross(z,x)
      R02.col(1) = (R02.col(1)) / (R02.col(1)).norm();
      R12=-R02;
    }
    if (Vc_before_1.norm()>1e-8) {
      R01.col(0) = Vc_before_1 / Vc_before_1.norm(); // 충돌직전 공의 지면에 대한 tangential상대속도의 방향을 x축으로 잡는다.
      R01.col(2) << 0, 0, 1; //ground 에서 contact의 z축은 world의 z축
      R01.col(1) = hat(R01.col(2)) * R01.col(0); // y축은 cross(z,x)
      R01.col(1) = (R01.col(1)) / (R01.col(1)).norm();
    }
    if (Vc_before_3.norm()>1e-8) {
      R13.col(0) = Vc_before_3 / Vc_before_3.norm(); // 충돌직전 공의 지면에 대한 tangential상대속도의 방향을 x축으로 잡는다.
      R13.col(2) << 0, 0, 1; //ground 에서 contact의 z축은 world의 z축
      R13.col(1) = hat(R13.col(2)) * R13.col(0); // y축은 cross(z,x)
      R13.col(1) = (R13.col(1)) / (R13.col(1)).norm();
    }

    ball[i].Jg.block(0, 0, 3, 3) = I_mat_3;
    ball[i].Jg.block(0, 3, 3, 3) = -hat(rc_01);
    ball[j].Jg.block(0, 0, 3, 3) = I_mat_3;
    ball[j].Jg.block(0, 3, 3, 3) = -hat(rc_13);
    ball[i].Jb.block(0, 0, 3, 3) = I_mat_3;
    ball[i].Jb.block(0, 3, 3, 3) = -hat(rc_02);
    ball[j].Jb.block(0, 0, 3, 3) = I_mat_3;
    ball[j].Jb.block(0, 3, 3, 3) = -hat(rc_12);

    bterm_0 << ball[i].mass * g, 0, 0, 0;
    bterm_0 += Ja.transpose() * (hat(ball[i].w) * (ball[i].M.block(3, 3, 3, 3) * ball[i].w));
    bterm_1 << ball[j].mass * g, 0, 0, 0;
    bterm_1 += Ja.transpose() * (hat(ball[j].w) * (ball[j].M.block(3, 3, 3, 3) * ball[j].w));

    Mapp = ((ball[i].Jb * (ball[i].M.inverse()) * ball[i].Jb.transpose()) +
            ((ball[j].Jb * (ball[j].M.inverse()) * ball[j].Jb.transpose()))).inverse();


    while (true) {
      std::cout<<"Vc_before1: "<<Vc_before_1<<std::endl;
      std::cout<<"Vc_before2: "<<Vc_before_2<<std::endl;
      std::cout<<"Vc_before3: "<<Vc_before_3<<std::endl;
      lambda1 << lx1, 0, lz1; //contact frame기준
      lambda2 << lx2, 0, lz2; //contact frame기준
      lambda3 << lx3, 0, lz3; //contact frame기준
      std::cout<<"lambda1: "<<lambda1<<std::endl;
      std::cout<<"lambda2: "<<lambda2<<std::endl;
      std::cout<<"lambda3: "<<lambda3<<std::endl;
      Vc_after_1= Vc_before_1 + (ball[0].Jg*(ball[0].M.inverse())*ball[0].Jg.transpose())*(R01*lambda1)-dt*(ball[0].Jg*ball[0].M.inverse()*bterm_0);
      Vc_after_2= Vc_before_2 + (ball[0].Jb*(ball[0].M.inverse())*ball[0].Jg.transpose())*(R01*lambda1)+Mapp.inverse()*R02*lambda2
                  -(ball[1].Jb*(ball[1].M.inverse())*ball[1].Jg.transpose())*(R13*lambda3)
                  -dt * (ball[i].Jb * ball[i].M.inverse() * bterm_0)
                  +dt * (ball[j].Jb * ball[j].M.inverse() * bterm_1);
      Vc_after_3= Vc_before_3 + (ball[1].Jg*(ball[1].M.inverse())*ball[1].Jg.transpose())*(ball[1].R*lambda1)-dt*(ball[1].Jg*ball[1].M.inverse()*bterm_1);
      prox_z1 = lz1 - alpha * (R01.transpose() * Vc_after_1)(2) * (1 / ball[0].M(2, 2)); //contact frame기준
      prox_z2 = lz2 - alpha * (R02.transpose() * Vc_after_2)(2) / (1 / Mapp(2, 2)); //contact frame기준
      prox_z3 = lz3 - alpha * (R13.transpose() * Vc_after_3)(2) * (1 / ball[1].M(2, 2)); //contact frame기준
      lz1= cons_prox_z(prox_z1);
      lz2= cons_prox_z(prox_z2);
      lz3= cons_prox_z(prox_z3);
      lambda1 << lx1, 0, lz1; //contact frame기준
      lambda2 << lx2, 0, lz2; //contact frame기준
      lambda3 << lx3, 0, lz3; //contact frame기준

      Vc_after_1= Vc_before_1 + (ball[0].Jg*(ball[0].M.inverse())*ball[0].Jg.transpose())*(R01*lambda1)-dt*(ball[0].Jg*ball[0].M.inverse()*bterm_0);
      Vc_after_2= Vc_before_2 + (ball[0].Jb*(ball[0].M.inverse())*ball[0].Jg.transpose())*(R01*lambda1)+Mapp.inverse()*R02*lambda2
                  -(ball[1].Jb*(ball[1].M.inverse())*ball[1].Jg.transpose())*(R13*lambda3)
                  -dt * (ball[i].Jb * ball[i].M.inverse() * bterm_0)
                  +dt * (ball[j].Jb * ball[j].M.inverse() * bterm_1);
      Vc_after_3= Vc_before_3 + (ball[1].Jg*(ball[1].M.inverse())*ball[1].Jg.transpose())*(ball[1].R*lambda1)-dt*(ball[1].Jg*ball[1].M.inverse()*bterm_1);

      prox_x1 = lx1 - beta * (R01.transpose() * Vc_after_1)(0) * (1 / ball[0].M(2, 2)); //contact frame기준
      prox_x2 = lx2 - beta * (R02.transpose() * Vc_after_2)(0) / (fmax(1 / Mapp(0,0),1/Mapp(1,1))); //contact frame기준
      prox_x3 = lx3 - beta * (R13.transpose() * Vc_after_3)(0) * (1 / ball[1].M(2, 2)); //contact frame기준

      lx1= cons_prox_x(prox_x1,lz1);
      lx2= cons_prox_x(prox_x2,lz2);
      lx3= cons_prox_x(prox_x3,lz3);

      lambda_prev<<lx_prev1,lx_prev2,lx_prev3,lz_prev1,lz_prev2,lz_prev3;
      lambda_now<<lx1,lx2,lx3,lz1,lz2,lz3;

      if ((lambda_prev-lambda_now).norm() < error_th) { //충분히 수렴한 경우
        lambda1 << lx1, 0, lz1; //contact frame기준
        lambda2 << lx2, 0, lz2; //contact frame기준
        lambda3 << lx3, 0, lz3; //contact frame기준

        ball[0].vel__ = ball[0].vel_ + (ball[0].M.inverse() * ((-bterm_0) * dt
                + ball[0].Jg.transpose() * R01 *lambda1+ball[0].Jb.transpose() * R02 *lambda2)).block(0, 0, 3,1);
        ball[0].w += (ball[0].M.inverse() * ((-bterm_0) * dt
                     + ball[0].Jg.transpose() * R01 *lambda1+ball[0].Jb.transpose() * R02 *lambda2)).block(3, 0, 3,1);
        ball[0].pos__ = ball[0].pos_ + dt * ball[0].vel__;

        ball[1].vel__ = ball[1].vel_ + (ball[1].M.inverse() * ((-bterm_1) * dt
                                                               + ball[1].Jg.transpose() * R13 *lambda3+ball[1].Jb.transpose() * R12 *lambda2)).block(0, 0, 3,1);
        ball[1].w += (ball[1].M.inverse() * ((-bterm_1) * dt
                                             + ball[1].Jg.transpose() * R13 *lambda3+ball[1].Jb.transpose() * R12 *lambda2)).block(3, 0, 3,1);
        ball[1].pos__ = ball[1].pos_ + dt * ball[1].vel__;

        break;
      }
      lx_prev1 = lx1;
      lz_prev1 = lz1;
      lx_prev2 = lx2;
      lz_prev2 = lz2;
      lx_prev3 = lx3;
      lz_prev3 = lz3;
    }
  }
  bool contact_ground(int i){
    if(ball[i].pos__(2) < ball[i].radius){
      return true;
    }
    return false;
  }
  bool contact_ball(int i, int j){
    Eigen::Vector3d posdiff=ball[i].pos__-ball[j].pos__;
    if(posdiff.norm()<(ball[i].radius+ball[j].radius)){
      return true;
    }
    return false;
  }
  int Contact_mode(int i){// return 0,1,2,3 for no contact, contact with ground, contact with ball, both contact case
      if(contact_ball(0,1)){ //ball끼리 만나는 경우
        if(contact_ground(i)){ // 공과 땅과 만나는 경우
          return 3;
        }
        else{
          return 2;
        }
      }
      if(contact_ground(i)){ //땅과만 만나는 경우
        return 1;
      }
      else { // contact이 없는 경우
        return 0;
      }

  }

  void freefall(int i){
      ball[i].vel__=ball[i].vel_+ball[i].acc*dt;
      ball[i].pos__=ball[i].pos_+ball[i].vel__*dt;
  }
  void setPosition(raisim::Visuals * sphere1, raisim::Visuals * sphere2) {
    sphere1->setPosition(ball[0].pos__);
    sphere2->setPosition(ball[1].pos__);
  }
 private:
  struct ball
  {
      double radius;
      double mass;
      Eigen::Vector3d pos_; //before impact
      Eigen::Vector3d pos__;//after impact
      Eigen::Vector3d vel_; //before impact
      Eigen::Vector3d vel__; //after impact
      Eigen::Vector3d acc;
      Eigen::Vector3d w;
      Eigen::MatrixXd Jg;//contact with ground 3*6
      Eigen::MatrixXd Jb;//contact with ball 3*6
      Eigen::Matrix3d R;//rotation matrix of body
      Eigen::MatrixXd M;//Mass matrix 6*6
      Eigen::Vector3d lambda_ground;
      Eigen::Vector3d lambda_ball;
      Eigen::Vector3d contact_point_b;
  };
  std::vector<ball> ball;
  double mu=0.8;
  double dt = 0.001;

  Eigen::Vector3d g={0,0,9.81};//TODO CHECK
  Eigen::MatrixXd Ja;
  Eigen::MatrixXd I_mat_6;
  Eigen::Matrix3d I_mat_3= ((Eigen::Matrix3d()<< Eigen::Matrix3d().setIdentity()).finished()); //identity matrix

  /// add state variables here
};

