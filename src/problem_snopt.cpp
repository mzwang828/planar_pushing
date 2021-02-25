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

  Problem nlp;
  Eigen::VectorXd initValues(n_step * 4);
  initValues.setZero();
  Eigen::VectorXd stateNominal(n_step*4+4), controlNominal((n_step-1)*3);
  for (int i = 0; i < n_step; ++i){
    // stateNominal.segment(i*4, 4) << 0.05*i*tStep, 0, 0, 0;
    stateNominal.segment(i*4, 4) << 0.15 * sin(0.2*i*tStep), 0.15 - 0.15*cos(0.2*i*tStep), 0.2*i*tStep, 0;
    initValues.segment(i*4, 4) << 0.15 * sin(0.2*i*tStep), 0.15 - 0.15*cos(0.2*i*tStep), 0.2*i*tStep, 0;
  }
  controlNominal.setZero();
  stateNominal.tail(4) << 0.15 * sin(0.2*n_step*tStep), 0.15 - 0.15*cos(0.2*n_step*tStep), 0.2*n_step*tStep, 0;

  initValues.head(4) << 0.0, 0.0, 0, 0;

  nlp.AddVariableSet(std::make_shared<ExVariables>(4*n_step, "state", initValues));
  nlp.AddVariableSet(std::make_shared<ExVariables>(3*(n_step-1), "control", controlNominal));
  nlp.AddConstraintSet(std::make_shared<ExConstraint>(10*(n_step-1)));
  nlp.AddCostSet(std::make_shared<ExCost>("cost", stateNominal.segment(0, n_step*4), controlNominal));

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

  getchar();

  variables.head(4) << 0.02,0,0,0;
  Problem nlp2;
  nlp2.AddVariableSet(std::make_shared<ExVariables>(4*n_step, "state", variables.segment(0, 4 * n_step)));
  nlp2.AddVariableSet(std::make_shared<ExVariables>(3*(n_step-1), "control", variables.segment(4 * n_step, 3 * (n_step-1))));
  nlp2.AddConstraintSet(std::make_shared<ExConstraint>(10*(n_step-1)));
  nlp2.AddCostSet(std::make_shared<ExCost>("cost", stateNominal.segment(4, 4*n_step), controlNominal));

  solver.Solve(nlp2);

  nlp2.PrintCurrent();

  
}