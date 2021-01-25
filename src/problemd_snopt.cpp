#include <iostream>
#include <fstream>
#include <ifopt/problem.h>
#include <ifopt/snopt_solver.h>
#include <ifopt/ipopt_solver.h>
#include "pusher/problemd.h"
#include <yaml-cpp/yaml.h>

using namespace ifopt;

int main()
{
  YAML::Node params = YAML::LoadFile("/home/mzwang/qsp_ws/src/pusher/Config/params.yaml");
  int n_step = params["n_step"].as<int>();         // number of steps

  Problem nlp;
  Eigen::VectorXd initStates(n_step * 4);
  initStates.setZero();

  Eigen::VectorXd stateNominal(n_step*4), controlNominal(n_step*2);
  for (int i = 0; i < n_step; ++i){
    stateNominal.segment(i*4, 4) << 0.05*i*0.03, 0,0,0;
    initStates.segment(i*4, 4) << 0.05*i*0.03, 0,0,0;
    controlNominal.segment(i*2, 2) << 0.05, 0;
  }
  initStates.head(4) << -0.0, 0.01, 15*3.14/180, 0;

  nlp.AddVariableSet(std::make_shared<ExVariables>(4*n_step, "state", initStates));
  nlp.AddVariableSet(std::make_shared<ExVariables>(4*n_step, "stateDot2"));
  nlp.AddVariableSet(std::make_shared<ExVariables>(4*n_step, "stateDot3"));
  nlp.AddVariableSet(std::make_shared<ExVariables>(2*n_step, "control"));
  nlp.AddConstraintSet(std::make_shared<ExConstraint>(4*5*(n_step-1)));

  nlp.AddCostSet(std::make_shared<ExCost>("cost", stateNominal, controlNominal));

  nlp.PrintCurrent();

  SnoptSolver solver;
  solver.Solve(nlp);

  // IpoptSolver solver;  
  // solver.SetOption("linear_solver", "mumps");
  // solver.SetOption("jacobian_approximation", "finite-difference-values");
  // solver.SetOption("max_cpu_time", 1e6);
  // solver.SetOption("max_iter", 30000);  
  // solver.Solve(nlp);

  nlp.PrintCurrent();

  Eigen::VectorXd variables = nlp.GetOptVariables()->GetValues();
  Eigen::Map<Eigen::MatrixXd> state(variables.segment(0, 4 * n_step).data(), 4, n_step);
  Eigen::Map<Eigen::MatrixXd> stateDot2(variables.segment(4 * 1 * n_step, 4 * n_step).data(), 4, n_step);
  Eigen::Map<Eigen::MatrixXd> stateDot3(variables.segment(4 * 2 * n_step, 4 * n_step).data(), 4, n_step);
  Eigen::Map<Eigen::MatrixXd> control(variables.segment(4 * 3 * n_step, 2 * n_step).data(), 2, n_step);

  std::cout.precision(5);

  std::cout << "state: \n" << state << "\n";
  std::cout << "control: \n" << control << "\n";
  std::cout << "stateDot2: \n" << stateDot2 << "\n";
  std::cout << "stateDot3: \n" << stateDot3 << "\n";

}