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

  Problem nlp;
  Eigen::VectorXd initValues(n_step * 4);
  initValues.setZero();
  Eigen::VectorXd stateNominal(n_step*4), controlNominal(n_step*3);
  for (int i = 0; i < n_step; ++i){
    stateNominal.segment(i*4, 4) << 0.05*i*0.03, 0,0,0;
    initValues.segment(i*4, 4) << 0.05*i*0.03, 0,0,0;
    controlNominal.segment(i*3, 3) << 0.0, 0, 0;
  }

  initValues.head(4) << 0.0, 0.1, 0.5, 0.0;

  nlp.AddVariableSet(std::make_shared<ExVariables>(4*n_step, "state", initValues));
  nlp.AddVariableSet(std::make_shared<ExVariables>(3*n_step, "control"));
  nlp.AddConstraintSet(std::make_shared<ExConstraint>(8*(n_step-1)));
  nlp.AddCostSet(std::make_shared<ExCost>("cost", stateNominal, controlNominal));

  nlp.PrintCurrent();

  SnoptSolver solver;
  solver.Solve(nlp);

  nlp.PrintCurrent();

  Eigen::VectorXd variables = nlp.GetOptVariables()->GetValues();
  Eigen::Map<Eigen::MatrixXd> state(variables.segment(0, 4 * n_step).data(), 4, n_step);
  Eigen::Map<Eigen::MatrixXd> control(variables.segment(4 * n_step, 3 * n_step).data(), 3, n_step);

  std::cout.precision(5);

  std::cout << "state: \n" << state << "\n";
  std::cout << "control: \n" << control << "\n";
  
}