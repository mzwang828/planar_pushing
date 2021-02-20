#include <ifopt/constraint_set.h>
#include <ifopt/cost_term.h>
#include <ifopt/variable_set.h>
#include <yaml-cpp/yaml.h>
#include <iostream>

typedef Eigen::Triplet<double> T;

namespace ifopt
{

  class ExVariables : public VariableSet
  {
  public:
    // Every variable set has a name, here "var_set1". this allows the constraints
    // and costs to define values and Jacobians specifically w.r.t this variable
    // set.
    ExVariables(int n) : ExVariables(n, "var_set1"){};
    ExVariables(int n, const std::string &name) : VariableSet(n, name)
    {
      xvar = Eigen::VectorXd::Zero(n);
      if (name == "state")
      {
        for (int i = 0; i < n; i++)
          xvar(i) = 0;
      }
      else if (name == "control")
      {
        for (int i = 0; i < n; i++)
        {
          xvar(i) = 0.05;
        }
      }
    }

    ExVariables(int n, const std::string &name, const Eigen::VectorXd &initValues) : VariableSet(n, name) 
    {
      xvar = Eigen::VectorXd::Zero(n);
      xvar = initValues;
      if (name == "state"){
        x0.push_back(initValues[0]);
        x0.push_back(initValues[1]);
        x0.push_back(initValues[2]);
        x0.push_back(initValues[3]);
      }
    }

    // Here is where you can transform the Eigen::Vector into whatever
    // internal representation of your variables you have (here two doubles, but
    // can also be complex classes such as splines, etc..
    void SetVariables(const VectorXd &x) override
    {
      for (int i = 0; i < x.size(); i++)
        xvar(i) = x(i);
    };

    // Here is the reverse transformation from the internal representation to
    // to the Eigen::Vector
    VectorXd GetValues() const override { return xvar; };

    // Each variable has an upper and lower bound set here
    VecBound GetBounds() const override
    {
      VecBound bounds(GetRows());
      if (GetName() == "state")
      {
        for (int i = 0; i < GetRows()/4; i++)
        {
          bounds.at(i*4) = Bounds(-10.0, 10.0);
          bounds.at(i*4 + 1) = Bounds(-1.0, 1.0);
          bounds.at(i*4 + 2) = Bounds(-3.14, 3.14);
          bounds.at(i*4 + 3) = Bounds(-0.35, 0.35);
        }
        bounds.at(0) = Bounds(x0[0], x0[0]);
        bounds.at(1) = Bounds(x0[1], x0[1]);
        bounds.at(2) = Bounds(x0[2], x0[2]);
        bounds.at(3) = Bounds(x0[3], x0[3]);
      }
      else if (GetName() == "control")
      {
        for (int i = 0; i < GetRows() / 3; i++)
        {
          bounds.at(3 * i) = Bounds(0.0, inf);
          bounds.at(3 * i + 1) = Bounds(-inf, inf);
          bounds.at(3 * i + 2) = Bounds(-inf, inf);
        }
      }
      return bounds;
    }

  private:
    Eigen::VectorXd xvar;
    std::vector<double> x0;
  };

  class ExConstraint : public ConstraintSet
  {
  public:
    ExConstraint(int n) : ExConstraint(n, "constraint1") {}

    // This constraint set just contains 1 constraint, however generally
    // each set can contain multiple related constraints.
    ExConstraint(int n, const std::string &name) : ConstraintSet(n, name)
    {
      YAML::Node params = YAML::LoadFile("/home/pengchang/build_ws/src/planar_pushing/Config/params.yaml");
      mu = params["mu"].as<double>();
      muGround = params["mu_g"].as<double>();
      m = params["m"].as<double>();
      xC = -params["length"].as<double>() / 2.0;
      mMax = params["mMax"].as<double>();
      fMax = muGround * m * 9.81;
      c = fMax / mMax;
      nStep = params["n_step"].as<int>();  
      tStep = params["t_step"].as<double>();
      
      L << 2/(fMax*fMax), 0, 0,
            0, 2/(fMax*fMax), 0,
            0, 0, 2/(mMax*mMax);
      
    }

    // The constraint value minus the constant value "1", moved to bounds.
    VectorXd GetValues() const override
    {
      VectorXd g(GetRows());
      VectorXd state = GetVariables()->GetComponent("state")->GetValues();
      VectorXd control = GetVariables()->GetComponent("control")->GetValues();
      for (int i = 1; i < nStep; i++)
      {
        double theta = state(i*4 + 2);
        double phi = state(i*4 + 3);
        double yC = - xC * tan(phi);
        double pdot = control(3*(i-1)+2);

        Eigen::Matrix3d R;
        R << cos(theta), -sin(theta), 0,
             sin(theta), cos(theta), 0,
             0, 0, 1;

        Eigen::MatrixXd J(2,3), B(3,2);
        J << 1, 0, -yC,
             0, 1, xC;

        B.col(0) = J.transpose() * Eigen::Vector2d(1, 0);
        B.col(1) = J.transpose() * Eigen::Vector2d(0, 1);


        Eigen::MatrixXd dynamics(4, 3), ddynamicsdTheta(4, 3), ddynamicsdPhi(4, 3);;
        dynamics.setZero();
        dynamics.topLeftCorner(3,2) = R*L*B;
        dynamics(3,2) = 1;


        g.segment(4 * (i - 1), 4) = (state.segment(4 * i, 4) - state.segment(4 * (i - 1), 4) -
                                     dynamics * control.segment(3*(i-1), 3) * tStep);

        // up
        g(4 * (nStep - 1) + i - 1) = -std::min(-pdot, 0.0) * (control(3 * (i-1) + 1) - mu * control(3 * (i-1)));

        // down
        g(5 * (nStep - 1) + i - 1) = -std::min(pdot, 0.0) * (control(3 * (i-1) + 1) + mu * control(3 * (i-1)));
        // friction cone
        g(6 * (nStep - 1) + i - 1) = control(3 * (i-1) + 1) - mu * control(3 * (i-1));
        g(7 * (nStep - 1) + i - 1) = control(3 * (i-1) + 1) + mu * control(3 * (i-1));

        // constraints for the pusher vel, v_n and v_t
        Eigen::MatrixXd Gc(2,3), dGcdPhi(2,3);;
        Gc.leftCols(2) = J * L * B;
        Gc.col(2) << 0.0, -xC/(cos(phi)*cos(phi));
        Eigen::Vector2d vPusher = Gc * control.segment(3*(i-1), 3);
        g(8 * (nStep - 1) + i - 1) = vPusher(0);   // normal velocity
        g(9 * (nStep - 1) + i - 1) = vPusher(1);   // tangential velocity
      }
      return g;
    };

    // The only constraint in this set is an equality constraint to 1.
    // Constant values should always be put into GetBounds(), not GetValues().
    // For inequality constraints (<,>), use Bounds(x, inf) or Bounds(-inf, x).
    VecBound GetBounds() const override
    {
      VecBound bounds(GetRows());
      for (int i = 0; i < 6 * (nStep - 1); ++i){
        bounds.at(i) = Bounds(0.0, 0.0);
      }
      for (int i = 0; i < nStep - 1; ++i){
        bounds.at(6*(nStep - 1)+i) = Bounds(-inf, 0.0);
        bounds.at(7*(nStep - 1)+i) = Bounds(0.0, inf);
        bounds.at(8*(nStep - 1)+i) = Bounds(0.0, 0.3);
        bounds.at(9*(nStep - 1)+i) = Bounds(-0.3, 0.3);
      }
      return bounds;
    }

    // This function provides the first derivative of the constraints.
    // In case this is too difficult to write, you can also tell the solvers to
    // approximate the derivatives by finite differences and not overwrite this
    // function, e.g. in ipopt.cc::use_jacobian_approximation_ = true
    void FillJacobianBlock(std::string var_set,
                           Jacobian &jac_block) const override
    {
      VectorXd state = GetVariables()->GetComponent("state")->GetValues();
      VectorXd control = GetVariables()->GetComponent("control")->GetValues();

      std::vector<T> tripletState, tripletControl;
      for (int i = 1; i < nStep; ++i)
      {            
        double theta = state(i * 4 + 2);
        double phi = state(i * 4 + 3);
        double yC = -xC * tan(phi);
        double pdot = control(3 * (i-1) + 2);

        Eigen::Matrix3d R;
        R << cos(theta), -sin(theta), 0,
            sin(theta), cos(theta), 0,
            0, 0, 1;

        Eigen::MatrixXd J(2, 3), B(3, 2);
        J << 1, 0, -yC,
            0, 1, xC;

        // B.col(0) = J.transpose() * Eigen::Vector2d(1, 0);
        // B.col(1) = J.transpose() * Eigen::Vector2d(0, 1);
        B = J.transpose();

        // Eigen::MatrixXd Gc(2,3), dGcdPhi(2,3);
        // Gc.leftCols(2) = J * L * B;
        // Gc.col(2) << 0.0, -xC/(cos(phi)*cos(phi));

        // Eigen::MatrixXd dRdTheta(3, 3), dBdPhi(3, 2), dJdPhi(2, 3);
        // dRdTheta << -sin(theta), -cos(theta), 0,
        //             cos(theta), -sin(theta), 0,
        //             0, 0, 0;

        // dBdPhi << 0, 0,
        //           0, 0,
        //           xC / (cos(phi) * cos(phi)), 0;

        // dJdPhi << 0, 0, xC / (cos(phi) * cos(phi)),
        //           0, 0, 0;

        // ddynamicsdTheta.setZero();
        // ddynamicsdTheta.topLeftCorner(3, 2) = dRdTheta * L * B;
        // ddynamicsdPhi.setZero();
        // ddynamicsdPhi.topLeftCorner(3, 2) = R * L * dBdPhi;

        // dGcdPhi.setZero();
        // dGcdPhi.leftCols(2) = dJdPhi * L * B + J * L * dBdPhi;
        // dGcdPhi.col(2) << 0.0, -2*xC*sin(phi)/(cos(phi)*cos(phi)*cos(phi));

        if (var_set == "state")
        {
          Eigen::MatrixXd ddynamicsdTheta(4, 3), ddynamicsdPhi(4, 3);
          Eigen::MatrixXd dGcdPhi(2,3);
          Eigen::MatrixXd dRdTheta(3, 3), dBdPhi(3, 2), dJdPhi(2, 3);

          dRdTheta << -sin(theta), -cos(theta), 0,
                      cos(theta), -sin(theta), 0,
                      0, 0, 0;

          dBdPhi << 0, 0,
                    0, 0,
                    xC / (cos(phi) * cos(phi)), 0;

          dJdPhi << 0, 0, xC / (cos(phi) * cos(phi)),
                    0, 0, 0;

          ddynamicsdTheta.setZero();
          ddynamicsdTheta.topLeftCorner(3, 2) = dRdTheta * L * B;
          ddynamicsdPhi.setZero();
          ddynamicsdPhi.topLeftCorner(3, 2) = R * L * dBdPhi;
          dGcdPhi.setZero();
          dGcdPhi.leftCols(2) = dJdPhi * L * B + J * L * dBdPhi;
          dGcdPhi.col(2) << 0.0, -2*xC*sin(phi)/(cos(phi)*cos(phi)*cos(phi));
          
          for (int j = 0; j < 4; ++j)
          {
            tripletState.push_back(T(4 * (i - 1) + j, 4 * i + j, 1));
            tripletState.push_back(T(4 * (i - 1) + j, 4 * (i - 1) + j, -1));
            // regarding theta
            tripletState.push_back(T(4 * (i - 1) + j, 4 * i + 2, (-ddynamicsdTheta * control.segment(3 * (i-1), 3) * tStep)(j)));
            // regarding phi
            tripletState.push_back(T(4 * (i - 1) + j, 4 * i + 3, (-ddynamicsdPhi * control.segment(3 * (i-1), 3) * tStep)(j)));
          }
          tripletState.push_back(T(8 * (nStep - 1) + i - 1, 4 * i + 3, (dGcdPhi * control.segment(3 * (i-1), 3))(0)));
          tripletState.push_back(T(9 * (nStep - 1) + i - 1, 4 * i + 3, (dGcdPhi * control.segment(3 * (i-1), 3))(1)));
        }
        else if (var_set == "control")
        {
          Eigen::MatrixXd dynamics(4, 3);
          dynamics.setZero();
          dynamics.topLeftCorner(3, 2) = R * L * B;
          dynamics(3, 2) = 1;

          Eigen::MatrixXd Gc(2,3);
          Gc.leftCols(2) = J * L * B;
          Gc.col(2) << 0.0, -xC/(cos(phi)*cos(phi));

          for (int j = 0; j < 4; ++j)
          {
            tripletControl.push_back(T(4 * (i - 1) + j, 3 * (i-1), -dynamics(j, 0) * tStep));
            tripletControl.push_back(T(4 * (i - 1) + j, 3 * (i-1) + 1, -dynamics(j, 1) * tStep));
            tripletControl.push_back(T(4 * (i - 1) + j, 3 * (i-1) + 2, -dynamics(j, 2) * tStep));
          }
          // up
          tripletControl.push_back(T(4 * (nStep - 1) + i - 1, 3 * (i-1), std::min(-pdot, 0.0) * mu));
          tripletControl.push_back(T(4 * (nStep - 1) + i - 1, 3 * (i-1) + 1, -std::min(-pdot, 0.0)));
          tripletControl.push_back(T(4 * (nStep - 1) + i - 1, 3 * (i-1) + 2, -pdot < 0.0 ? (control(3 * (i-1) + 1) - mu * control(3 * (i-1))) : 0));
          // down
          tripletControl.push_back(T(5 * (nStep - 1) + i - 1, 3 * (i-1), -std::min(pdot, 0.0) * mu));
          tripletControl.push_back(T(5 * (nStep - 1) + i - 1, 3 * (i-1) + 1, -std::min(pdot, 0.0)));
          tripletControl.push_back(T(5 * (nStep - 1) + i - 1, 3 * (i-1) + 2, pdot < 0.0 ? -(control(3 * (i-1) + 1) + mu * control(3 * (i-1))) : 0));
          // friction cone
          tripletControl.push_back(T(6 * (nStep - 1) + i - 1, 3 * (i-1), -mu));
          tripletControl.push_back(T(6 * (nStep - 1) + i - 1, 3 * (i-1) + 1, 1));
          tripletControl.push_back(T(7 * (nStep - 1) + i - 1, 3 * (i-1), mu));
          tripletControl.push_back(T(7 * (nStep - 1) + i - 1, 3 * (i-1) + 1, 1));
          // pusher velocity
          tripletControl.push_back(T(8 * (nStep - 1) + i - 1, 3 * (i-1), Gc(0,0)));
          tripletControl.push_back(T(8 * (nStep - 1) + i - 1, 3 * (i-1)+1, Gc(0,1)));
          tripletControl.push_back(T(8 * (nStep - 1) + i - 1, 3 * (i-1)+2, Gc(0,2)));
          tripletControl.push_back(T(9 * (nStep - 1) + i - 1, 3 * (i-1), Gc(1,0)));
          tripletControl.push_back(T(9 * (nStep - 1) + i - 1, 3 * (i-1)+1, Gc(1,1)));
          tripletControl.push_back(T(9 * (nStep - 1) + i - 1, 3 * (i-1)+2, Gc(1,2)));
        }
      }
      if (var_set == "state")
      {
        jac_block.setFromTriplets(tripletState.begin(), tripletState.end());
      }
      else if (var_set == "control")
      {
        jac_block.setFromTriplets(tripletControl.begin(), tripletControl.end());
      }
    }

  private:
    double xC, fMax, mMax, c, mu, m, muGround, tStep;

    Eigen::Matrix3d L;

    int nStep;
  };

  class ExCost : public CostTerm
  {
  public:
    ExCost() : ExCost("cost_term1") {}
    ExCost(const std::string &name) : CostTerm(name) {}
    ExCost(const std::string &name,
           const Eigen::VectorXd &stateIn,
           const Eigen::VectorXd &controlIn) : CostTerm(name)
    {
      stateNominal = stateIn;
      controlNominal = controlIn;
      YAML::Node params = YAML::LoadFile("/home/pengchang/build_ws/src/planar_pushing/Config/params.yaml");
      QFinal = params["QFinal"].as<double>();
      Q = params["Q"].as<double>();
      R = params["R"].as<double>();
      nStep = params["n_step"].as<int>();  
      weightQ.resize(4, 4);
      weightR.resize(3, 3);
      weightQ = Eigen::Vector4d(3, 3, 0.1, 0).asDiagonal();
      weightR = Eigen::Vector3d(1, 1, 0.01).asDiagonal();
    }

    double GetCost() const override
    {
      VectorXd state = GetVariables()->GetComponent("state")->GetValues();
      VectorXd control = GetVariables()->GetComponent("control")->GetValues();

      double cost = 0;
      for (int i = 0; i < nStep - 1; ++i)
      {
        cost += (state.segment(i * 4, 4) - stateNominal.segment(i * 4, 4)).transpose() *
                Q * weightQ * (state.segment(i * 4, 4) - stateNominal.segment(i * 4, 4));
        cost += (control.segment(i * 3, 3) - controlNominal.segment(i * 3, 3)).transpose() *
                R * weightR * (control.segment(i * 3, 3) - controlNominal.segment(i * 3, 3));
      }
      cost += (state.tail(4) - stateNominal.tail(4)).transpose() *
              QFinal * weightQ * (state.tail(4) - stateNominal.tail(4));
      return cost;
    }

    void FillJacobianBlock(std::string var_set, Jacobian &jac) const override
    {
      VectorXd state = GetVariables()->GetComponent("state")->GetValues();
      VectorXd control = GetVariables()->GetComponent("control")->GetValues();
      std::vector<T> tripletState, tripletControl;

      if (var_set == "state")
      {
        for (int i = 0; i < nStep - 1; ++i)
        {
          tripletState.push_back(T(0, i * 4, 2 * Q * weightQ(0, 0) * (state(i * 4) - stateNominal(i * 4))));
          tripletState.push_back(T(0, i * 4 + 1, 2 * Q * weightQ(1, 1) * (state(i * 4 + 1) - stateNominal(i * 4 + 1))));
          tripletState.push_back(T(0, i * 4 + 2, 2 * Q * weightQ(2, 2) * (state(i * 4 + 2) - stateNominal(i * 4 + 2))));
          tripletState.push_back(T(0, i * 4 + 3, 2 * Q * weightQ(3, 3) * (state(i * 4 + 3) - stateNominal(i * 4 + 3))));
        }
        tripletState.push_back(T(0, (nStep - 1) * 4, 2 * QFinal * weightQ(0, 0) * (state((nStep - 1) * 4) - stateNominal((nStep - 1) * 4))));
        tripletState.push_back(T(0, (nStep - 1) * 4 + 1, 2 * QFinal * weightQ(1, 1) * (state((nStep - 1) * 4 + 1) - stateNominal((nStep - 1) * 4 + 1))));
        tripletState.push_back(T(0, (nStep - 1) * 4 + 2, 2 * QFinal * weightQ(2, 2) * (state((nStep - 1) * 4 + 2) - stateNominal((nStep - 1) * 4 + 2))));
        tripletState.push_back(T(0, (nStep - 1) * 4 + 3, 2 * QFinal * weightQ(3, 3) * (state((nStep - 1) * 4 + 3) - stateNominal((nStep - 1) * 4 + 3))));
        jac.setFromTriplets(tripletState.begin(), tripletState.end());
      }
      else if (var_set == "control")
      {
        for (int i = 0; i < nStep - 1; ++i)
        {
          tripletControl.push_back(T(0, i * 3, 2 * R * weightR(0, 0) * (control(i * 3) - controlNominal(i * 3))));
          tripletControl.push_back(T(0, i * 3 + 1, 2 * R * weightR(1, 1) * (control(i * 3 + 1) - controlNominal(i * 3 + 1))));
          tripletControl.push_back(T(0, i * 3 + 2, 2 * R * weightR(2, 2) * (control(i * 3 + 2) - controlNominal(i * 3 + 2))));
        }
        jac.setFromTriplets(tripletControl.begin(), tripletControl.end());
      }
    }

  private:
    Eigen::VectorXd stateNominal, controlNominal;
    Eigen::MatrixXd weightQ, weightR;
    double QFinal, Q, R;
    int nStep;
  };

} // namespace ifopt
