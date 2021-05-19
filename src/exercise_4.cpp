#include "raisim/RaisimServer.hpp"
#include "exercise_4_STUDENTID.hpp"


int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);

  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world); // visualization server
  world.addGround();

  // kinova
  auto kinova = world.addArticulatedSystem(binaryPath.getDirectory() + "/rsc/kinova/urdf/kinova_no_fingers.urdf");
  kinova->setName("kinova");
  server.focusOn(kinova);

  // kinova configuration
  Eigen::VectorXd gc(kinova->getGeneralizedCoordinateDim()), gv(kinova->getDOF());
  gc << 0.0, 2.56, -1.5, 0.0, 2.0, 0.0;
  gv << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6;
  kinova->setState(gc, gv);

  /// this function updates internal variables in raisim (such as the mass matrix and nonlinearities)
  world.integrate1();

  /// DO NOT USE OR MODIFY THIS CODE. /////////////////////////////////////////////

  /// You can use these functions to debug but NOT IN YOUR SUBMISSION
//   std::cout << static_cast<Eigen::Vector3d>(kinova->getBodyCOM_W()) << std::endl;
//   std::cout << kinova->getDenseJacobian() << std::endl;
//    auto debugSphere = server.addVisualSphere("debug_sphere", 0.15);
//    debugSphere->setColor(1,0,0,1);
//    debugSphere->setPosition(getEndEffectorPosition(gc));
//  std::cout << "driven Mass Matrix : " << getMassMatrix(gc) << "\n" << std::endl;
//  std::cout << "raisim Mass Matrix : " << kinova->getMassMatrix() << "\n" << std::endl;
//  std::cout << (getMassMatrix(gc) - kinova->getMassMatrix().e()).norm() << "\n" << std::endl;
  if((getMassMatrix(gc) - kinova->getMassMatrix().e()).norm() < 1e-8)
    std::cout<<"passed "<<std::endl;

//  std::cout << "driven C : " << getNonlinearities(gc, gv) << "\n" << std::endl;
//  std::cout << "raisim C : " << kinova->getNonlinearities() << "\n" << std::endl;
//  std::cout << (getNonlinearities(gc, gv) - kinova->getNonlinearities().e()).norm() << std::endl;
  if((getNonlinearities(gc, gv) - kinova->getNonlinearities().e()).norm() < 1e-8)
    std::cout<<"passed "<<std::endl;

  // visualization
  server.launchServer();

  /// visualize to debug!
  for (int i=0; i<2000000; i++)
    std::this_thread::sleep_for(std::chrono::microseconds(1000));

  server.killServer();
}
