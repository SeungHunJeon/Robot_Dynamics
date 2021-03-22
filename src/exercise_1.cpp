// This file is part of RaiSim. You must obtain a valid license from RaiSim Tech
// Inc. prior to usage.

#include "raisim/RaisimServer.hpp"
#include "exercise_1_STUDENTID.hpp"


int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);

  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world); // visualization server
  world.addGround();

  // kinova
  auto kinova = world.addArticulatedSystem(binaryPath.getDirectory() + "/rsc/kinova/urdf/kinova.urdf");
  kinova->setName("kinova");
  server.focusOn(kinova);

  // kinova configuration
  Eigen::VectorXd jointNominalConfig(kinova->getGeneralizedCoordinateDim());
  jointNominalConfig << 0.0, -1.547, 1,2, -0.5, 1.1, 0.0, 0.0, 0.0, 0.0;
  kinova->setGeneralizedCoordinate(jointNominalConfig);

  raisim::Vec<3> position1;
  kinova->getBodyPosition(4, position1);
  kinova->getDOF();
  std::cout << position1 << "\n";

  // debug sphere
  auto debugSphere = server.addVisualSphere("debug_sphere", 0.15);
  debugSphere->setColor(1,0,0,1);
  debugSphere->setPosition(getEndEffectorPosition(jointNominalConfig));

  // visualization
  server.launchServer(9000);
  for (int i=0; i<2000000; i++)
    std::this_thread::sleep_for(std::chrono::microseconds(1000));

  server.killServer();
}
