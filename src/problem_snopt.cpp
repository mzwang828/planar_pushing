#include <iostream>
#include <fstream>
#include <ifopt/problem.h>
#include <ifopt/snopt_solver.h>
#include "pusher/problem.h"
#include <yaml-cpp/yaml.h>

using namespace ifopt;

int main()
{
  YAML::Node params = YAML::LoadFile("/home/mzwang/qsp_ws/src/pusher/Config/params.yaml");
  int n_step = params["n_step"].as<int>();         // number of steps
  double tStep = params["t_step"].as<double>();


  // control for 8 shape
  // float radius = 0.15;
  // float targetAV = 0.2;
  // int nPointsPi = 3.14/(targetAV*tStep);      
  // Eigen::VectorXd eightNomi;  
  // eightNomi.resize(nPointsPi*4*4);
  // for (int i = 0; i < nPointsPi; ++i){
  //     eightNomi.segment(i*4, 4) << radius * sin(targetAV*i*tStep), radius - radius*cos(targetAV*i*tStep), targetAV*i*tStep, 0;
  //     eightNomi.segment(4*nPointsPi + i*4, 4) << -radius * sin(targetAV*i*tStep), 3*radius - radius*cos(targetAV*i*tStep), 3.14 - targetAV*i*tStep, 0;
  //     eightNomi.segment(2*4*nPointsPi + i*4, 4) << radius * sin(targetAV*i*tStep), 3*radius + radius*cos(targetAV*i*tStep), - targetAV*i*tStep, 0;
  //     eightNomi.segment(3*4*nPointsPi + i*4, 4) << -radius * sin(targetAV*i*tStep), radius + radius*cos(targetAV*i*tStep), -3.14 + targetAV*i*tStep, 0;
  // }
  // n_step = nPointsPi * 4;

  
  Problem nlp;
  Eigen::VectorXd initValues(n_step * 4);
  initValues.setZero();
  Eigen::VectorXd stateNominal(n_step*4), controlNominal((n_step-1)*3);
  for (int i = 0; i < n_step; ++i){
    stateNominal.segment(i*4, 4) << 0.05*i*tStep, 0, 0, 0;
    // stateNominal.segment(i*4, 4) << 0.15 * sin(0.2*i*tStep), 0.15 - 0.15*cos(0.2*i*tStep), 0.2*i*tStep, 0;
    // initValues.segment(i*4, 4) << 0.15 * sin(0.2*i*tStep), 0.15 - 0.15*cos(0.2*i*tStep), 0.2*i*tStep, 0;
  }
  controlNominal.setZero();
  for (int i = 0; i < n_step-1; ++i){
    controlNominal.segment(i*3, 3) << 0.305, 0, 0;
  }
  // stateNominal.tail(4) << 0.15 * sin(0.2*n_step*tStep), 0.15 - 0.15*cos(0.2*n_step*tStep), 0.2*n_step*tStep, 0;

  initValues.head(4) << 0, 0.03, 0, 0;
  nlp.AddVariableSet(std::make_shared<ExVariables>(4*n_step, "state", initValues));
  nlp.AddVariableSet(std::make_shared<ExVariables>(3*(n_step-1), "control", controlNominal));
  nlp.AddConstraintSet(std::make_shared<ExConstraint>(10*(n_step-1)));
  nlp.AddCostSet(std::make_shared<ExCost>("cost", stateNominal, controlNominal));

  nlp.PrintCurrent();

  SnoptSolver solver;
  solver.Solve(nlp);

  nlp.PrintCurrent();

  Eigen::VectorXd variables = nlp.GetOptVariables()->GetValues();
  Eigen::Map<Eigen::MatrixXd> state(variables.segment(0, 4 * n_step).data(), 4, n_step);
  Eigen::Map<Eigen::MatrixXd> control(variables.segment(4 * n_step, 3 * (n_step-1)).data(), 3, n_step-1);
  Eigen::Map<Eigen::MatrixXd> reference(stateNominal.data(), 4, n_step);


  std::cout.precision(5);
  std::cout << "reference: \n" << reference << "\n";
  std::cout << "state: \n" << state << "\n";
  std::cout << "control: \n" << control << "\n";

  // std::string trajPathStr = "/home/mzwang/qsp_ws/src/pusher/logs/state.txt";
  // const char* trajPath = trajPathStr.c_str();
  // std::ofstream trajFile;
  // trajFile.open(trajPath);
  // if (trajFile.is_open()){
  //   for(int i=0; i<n_step; i++){
  //     trajFile << state(0, i) << ",";        
  //   }
  //   trajFile << "\n";
  //   for(int i=0; i<n_step; i++){
  //     trajFile << state(1, i) << ",";        
  //   }
  // }
  // else{
  //   std::cout << " WARNING: Unable to open the UR trajectory file.\n";
  // }
  // trajFile.close();

}