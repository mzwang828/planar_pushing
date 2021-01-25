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
        for (int i = 0; i < n / 2; i++)
        {
          xvar(2 * i) = 0.05;
          xvar(2 * i + 1) = 0.0;
        }
      }
      else if (name == "stateDot2")
      {
        for (int i = 0; i < n; i++)
          xvar(i) = 0;
      }
      else if (name == "stateDot3")
      {
        for (int i = 0; i < n; i++)
          xvar(i) = 0;
      }
    }

    ExVariables(int n, const std::string &name, const Eigen::VectorXd &initValues) : ExVariables(n, name)
    {
      xvar = Eigen::VectorXd::Zero(n);
      xvar = initValues;
      if (name == "state")
      {
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
        for (int i = 0; i < GetRows() / 4; i++)
        {
          bounds.at(i * 4) = Bounds(-1.0, 1.0);
          bounds.at(i * 4 + 1) = Bounds(-1.0, 1.0);
          bounds.at(i * 4 + 2) = Bounds(-3.14, 3.14);
          bounds.at(i * 4 + 3) = Bounds(-0.05, 0.05);
        }
        bounds.at(0) = Bounds(x0[0], x0[0]);
        bounds.at(1) = Bounds(x0[1], x0[1]);
        bounds.at(2) = Bounds(x0[2], x0[2]);
        bounds.at(3) = Bounds(x0[3], x0[3]);
      }
      else if (GetName() == "control")
      {
        for (int i = 0; i < GetRows(); i++)
          bounds.at(i) = Bounds(-0.1, 0.1);
      }
      else if (GetName() == "stateDot2")
      {
        for (int i = 0; i < GetRows(); i++)
          bounds.at(i) = Bounds(-0.1, 0.1);
      }
      else if (GetName() == "stateDot3")
      {
        for (int i = 0; i < GetRows(); i++)
          bounds.at(i) = Bounds(-0.1, 0.1);
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
      YAML::Node params = YAML::LoadFile("/home/mzwang/qsp_ws/src/pusher/Config/params.yaml");
      mu = params["mu"].as<double>();
      muGround = params["mu_g"].as<double>();
      m = params["m"].as<double>();
      px = -params["length"].as<double>() / 2.0;
      mMax = params["mMax"].as<double>();
      fMax = muGround * m * 9.81;
      c = fMax / mMax;
      nStep = n / 20 + 1;
      tStep = params["t_step"].as<double>();
    }

    // The constraint value minus the constant value "1", moved to bounds.
    VectorXd GetValues() const override
    {
      VectorXd g(GetRows());
      VectorXd state = GetVariables()->GetComponent("state")->GetValues();
      VectorXd control = GetVariables()->GetComponent("control")->GetValues();
      VectorXd stateDot2 = GetVariables()->GetComponent("stateDot2")->GetValues();
      VectorXd stateDot3 = GetVariables()->GetComponent("stateDot3")->GetValues();

      for (int i = 1; i < nStep; i++)
      {
        // get the velocity at step n+1
        double py = state(i * 4 + 3);
        double yt = (mu * c * c - px * py + mu * px * px) / (c * c + py * py - mu * px * py);
        double yb = (-mu * c * c - px * py - mu * px * px) / (c * c + py * py + mu * px * py);
        Eigen::Matrix2d C, Q, P1, P2, P3;
        C << cos(state(i * 4 + 2)), sin(state(i * 4 + 2)),
            -sin(state(i * 4 + 2)), cos(state(i * 4 + 2));
        Q << c * c + px * px, px * py, px * py, c * c + py * py;
        double fenmu = (c * c + px * px + py * py);
        Q = Q / fenmu;
        P1.setIdentity();
        P2 << 0.0, 0.0, yt, -1.0;
        P3 << 0.0, 0.0, yb, -1.0;

        Eigen::MatrixXd dynamicsStick(4, 2), dynamicsUp(4, 2), dynamicsDown(4, 2);
        dynamicsStick.topRows(2) = C.transpose() * Q * P1;
        dynamicsStick.row(2) << -py / fenmu, px;
        dynamicsStick.row(3) << 0.0, 0.0;

        dynamicsUp.topRows(2) = C.transpose() * Q * P2;
        dynamicsUp.row(2) << (yt * px) / fenmu, -px;
        dynamicsUp.row(3) << -yt, 0.0;

        dynamicsDown.topRows(2) = C.transpose() * Q * P3;
        dynamicsDown.row(2) << (yb * px) / fenmu, -px;
        dynamicsDown.row(3) << -yb, 0.0;

        Eigen::VectorXd stateDot(4);
        stateDot = dynamicsStick * control.segment(2 * i, 2) + stateDot2.segment(4 * i, 4) + stateDot3.segment(4 * i, 4);

        g.segment(4 * (i - 1), 4) = state.segment(4 * i, 4) - state.segment(4 * (i - 1), 4) - stateDot * tStep;

        // std::cout << std::scientific << g.segment(4 * (i - 1), 4).transpose() << "\n";
        // up
        g.segment(4 * (nStep - 1) + 4 * (i - 1), 4) = -std::min(-control(i * 2 + 1) + yt * control(i * 2), 0.0) *100 *
                                                      (stateDot2.segment(4 * i, 4) - dynamicsUp * control.segment(2 * i, 2));

        g.segment(4 * 2 * (nStep - 1) + 4 * (i - 1), 4) = -std::min(control(i * 2 + 1) - yt * control(i * 2), 0.0) *100*
                                                          (stateDot2.segment(4 * i, 4));

        // down
        g.segment(4 * 3 * (nStep - 1) + 4 * (i - 1), 4) = -std::min(control(i * 2 + 1) - yb * control(i * 2), 0.0) * 100*
                                                          (stateDot3.segment(4 * i, 4) - dynamicsDown * control.segment(2 * i, 2));

        g.segment(4 * 4 * (nStep - 1) + 4 * (i - 1), 4) = -std::min(-control(i * 2 + 1) + yb * control(i * 2), 0.0) *100*
                                                          (stateDot3.segment(4 * i, 4));

        if (control(i * 2 + 1) > yt * control(i * 2)){
          std::cout << "UP, ";
          // std::cout << std::scientific <<(stateDot2.segment(4 * i, 4) - dynamicsUp * control.segment(2 * i, 2)).transpose() << "), ";
        } else if(control(i * 2 + 1) < yb * control(i * 2)){
          std::cout << "DOWN ( ";
        } else {
          std::cout << "STICK, ";
        }
      }
      std::cout << "\n";

      return g;
    };

    // The only constraint in this set is an equality constraint to 1.
    // Constant values should always be put into GetBounds(), not GetValues().
    // For inequality constraints (<,>), use Bounds(x, inf) or Bounds(-inf, x).
    VecBound GetBounds() const override
    {
      VecBound bounds(GetRows());
      for (int i = 0; i < GetRows(); i++)
        bounds.at(i) = Bounds(0.0, 0.0);
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
      VectorXd stateDot2 = GetVariables()->GetComponent("stateDot2")->GetValues();
      VectorXd stateDot3 = GetVariables()->GetComponent("stateDot3")->GetValues();

      std::vector<T> tripletState, tripletControl, tripletStateDot2, tripletStateDot3;
      for (int i = 1; i < nStep; ++i)
      {
        double py = state(i * 4 + 3);
        double yt = (mu * c * c - px * py + mu * px * px) / (c * c + py * py - mu * px * py);
        double yb = (-mu * c * c - px * py - mu * px * px) / (c * c + py * py + mu * px * py);
        Eigen::Matrix2d C, Q, dC, dQ, P1, P2, P3, dP2, dP3;
        C << cos(state(i * 4 + 2)), sin(state(i * 4 + 2)),
            -sin(state(i * 4 + 2)), cos(state(i * 4 + 2));
        Q << c * c + px * px, px * py, px * py, c * c + py * py;
        double fenmu = (c * c + px * px + py * py);
        Q = Q / fenmu;
        P1.setIdentity();
        P2 << 0.0, 0.0, yt, -1.0;
        P3 << 0.0, 0.0, yb, -1.0;
        // dCdtheta
        dC << -sin(state(i * 4 + 2)), cos(state(i * 4 + 2)),
             -cos(state(i * 4 + 2)), -sin(state(i * 4 + 2));
        // dQdpy
        dQ << 2*py*(c*c+px*px), px*(c*c+px*px-py*py), px*(c*c+px*px-py*py), 2*px*px*py;
        dQ = dQ / (fenmu*fenmu);
        dP2 << 0.0,0.0,(c*c*mu*mu*px-2*c*c*mu*py-c*c*px+mu*mu*px*px*px-2*mu*px*px*py+px*py*py)/((-c*c+mu*px*py-py*py) * (-c*c+mu*px*py-py*py)), 0.0;
        dP3 << 0.0,0.0,(c*c*mu*mu*px+2*c*c*mu*py-c*c*px+mu*mu*px*px*px+2*mu*px*px*py+px*py*py)/((c*c+mu*px*py+py*py) * (c*c+mu*px*py+py*py)), 0.0;

        Eigen::Vector2d dStickdtheta, dUpdtheta, dDowndtheta;
        dStickdtheta = dC.transpose() * Q * control.segment(2 * i, 2);
        dUpdtheta = dC.transpose()*Q*P2*control.segment(2 * i, 2);
        dDowndtheta = dC.transpose()*Q*P3*control.segment(2 * i, 2);

        Eigen::MatrixXd dynamicsStick(4, 2), dynamicsUp(4, 2), dynamicsDown(4, 2);
        dynamicsStick.topRows(2) = C.transpose() * Q * P1;
        dynamicsStick.row(2) << -py / fenmu, px;
        dynamicsStick.row(3) << 0.0, 0.0;

        dynamicsUp.topRows(2) = C.transpose() * Q * P2;
        dynamicsUp.row(2) << (yt * px) / fenmu, -px;
        dynamicsUp.row(3) << -yt, 0.0;

        dynamicsDown.topRows(2) = C.transpose() * Q * P3;
        dynamicsDown.row(2) << (yb * px) / fenmu, -px;
        dynamicsDown.row(3) << -yb, 0.0;

        Eigen::Vector4d dStickdpy, dUpdpy, dDowndpy;
        dStickdpy.head(2) = C.transpose() * dQ * control.segment(2 * i, 2);
        dStickdpy(2) = control(2*i)*(c*c + px*px - py*py)/(fenmu * fenmu);
        dStickdpy(3) = 0.0;

        dUpdpy.head(2) = (C.transpose()*dQ*P2 + C.transpose()*Q*dP2)*control.segment(2*i,2);
        dUpdpy(2) = control(2*i)*(-2*px*py*(c*c*mu+mu*px*px-px*py)/(fenmu*fenmu*(c*c-mu*px*py+py*py)) - px*(2*py-mu*px)*(c*c*mu+mu*px*px-px*py)/(fenmu*(c*c-mu*px*py+py*py)*(c*c-mu*px*py+py*py)) - px*px/(fenmu*(c*c-mu*px*py+py*py)));
        dUpdpy(3) = -control(2*i)*(dP2(1,0));

        dDowndpy.head(2) = (C.transpose()*dQ*P3 + C.transpose()*Q*dP3)*control.segment(2*i,2);
        dDowndpy(2) = control(2*i)*(-px*px/((fenmu)*(c*c+mu*px*py+py*py)) - 2*px*py*(-mu*c*c-mu*px*px-px*py)/(fenmu*fenmu*(c*c+mu*px*py+py*py))-px*(mu*px+2*py)*(-c*c*mu-mu*px*px-px*py)/(fenmu*(c*c+mu*px*py+py*py)*(c*c+mu*px*py+py*py)));
        dDowndpy(3) = -control(2*i)*(dP3(1,0));

        for (int j = 0; j < 4; ++j)
        {
          if (var_set == "state")
          {
            tripletState.push_back(T(4 * (i - 1) + j, 4 * i + j, 1));
            tripletState.push_back(T(4 * (i - 1) + j, 4 * (i-1) + j, -1));
            // df1/dpy
            tripletState.push_back(T(4 * (i - 1) + j, 4 * i + 3, - dStickdpy(j)*tStep));
            // df2/dpy
            tripletState.push_back(T(4 * (nStep - 1) + 4 * (i - 1) + j, 4*i+3, -control(i*2+1)+yt*control(i*2) < 0? -dP2(1,0)*control(i*2)*(stateDot2.segment(4 * i, 4) - dynamicsUp * control.segment(2 * i, 2))(j) -
                                                                                                                    (-control(i * 2 + 1) + yt * control(i * 2)) * (-dUpdpy(j))*100 :0.0));
            tripletState.push_back(T(4 * 2 * (nStep - 1) + 4 * (i - 1) + j, 4*i+3, control(i*2+1)+yt*control(i*2) < 0? dP2(1,0)*control(i*2)*stateDot2.segment(4 * i, 4)(j)*100 : 0.0));
            // df3/dpy
            tripletState.push_back(T(4 * 3 * (nStep - 1) + 4 * (i - 1) + j, 4*i+3, control(i * 2 + 1) - yb * control(i * 2) < 0? dP3(1,0)*control(i*2)*(stateDot3.segment(4 * i, 4) - dynamicsDown * control.segment(2 * i, 2))(j) -
                                                                                                                                 (control(i * 2 + 1) - yb * control(i * 2)) * (-dDowndpy(j))*100: 0.0));
            tripletState.push_back(T(4 * 4 * (nStep - 1) + 4 * (i - 1) + j, 4*i+3, -control(i * 2 + 1) + yb * control(i * 2) < 0? -dP3(1,0)*control(i*2)*stateDot3.segment(4 * i, 4)(j)*100 : 0.0));
          }
          else if (var_set == "control")
          {
            tripletControl.push_back(T(4 * (i - 1) + j, 2 * i, - dynamicsStick(j,0) * tStep));
            tripletControl.push_back(T(4 * (i - 1) + j, 2 * i + 1, - dynamicsStick(j,1) * tStep));
            
            tripletControl.push_back(T(4 * (nStep - 1) + 4 * (i - 1) + j, 2*i, -control(i*2+1)+yt*control(i*2) < 0? -yt*(stateDot2.segment(4 * i, 4) - dynamicsUp * control.segment(2 * i, 2))(j) -
                                                                                                                    (-control(i * 2 + 1) + yt * control(i * 2)) * (-dynamicsUp(j,0))*100:0.0));
            tripletControl.push_back(T(4 * (nStep - 1) + 4 * (i - 1) + j, 2*i+1, -control(i*2+1)+yt*control(i*2) < 0? (stateDot2.segment(4 * i, 4) - dynamicsUp * control.segment(2 * i, 2))(j) - 
                                                                                                                    (-control(i * 2 + 1) + yt * control(i * 2)) * dynamicsUp(j,1)*100:0.0));
            tripletControl.push_back(T(4 * 2 * (nStep - 1) + 4 * (i - 1) + j, 2*i, control(i*2+1)+yt*control(i*2) < 0? yt*stateDot2(4 * i + j)*100 : 0.0));
            tripletControl.push_back(T(4 * 2 * (nStep - 1) + 4 * (i - 1) + j, 2*i+1, control(i*2+1)+yt*control(i*2) < 0? -stateDot2(4 * i + j)*100 : 0.0));

            tripletControl.push_back(T(4 * 3 * (nStep - 1) + 4 * (i - 1) + j, 2*i, control(i * 2 + 1) - yb * control(i * 2) < 0? yb*(stateDot3.segment(4 * i, 4) - dynamicsDown * control.segment(2 * i, 2))(j) -
                                                                                                                                 (control(i * 2 + 1) - yb * control(i * 2)) * (-dynamicsDown(j,0))*100 : 0.0));
            tripletControl.push_back(T(4 * 3 * (nStep - 1) + 4 * (i - 1) + j, 2*i+1, control(i * 2 + 1) - yb * control(i * 2) < 0? -(stateDot3.segment(4 * i, 4) - dynamicsDown * control.segment(2 * i, 2))(j) - 
                                                                                                                                    (control(i * 2 + 1) - yb * control(i * 2)) * dynamicsDown(j,1)*100:0.0));
            tripletControl.push_back(T(4 * 4 * (nStep - 1) + 4 * (i - 1) + j, 2*i, -control(i * 2 + 1) + yb * control(i * 2) < 0? -yb*stateDot3(4+i+j)*100 : 0.0));
            tripletControl.push_back(T(4 * 4 * (nStep - 1) + 4 * (i - 1) + j, 2*i+1, -control(i * 2 + 1) + yb * control(i * 2) < 0? stateDot3(4+i+j)*100 : 0.0));
          }
          else if (var_set == "stateDot2")
          {
            tripletStateDot2.push_back(T(4 * (i - 1) + j, 4 * i + j, -tStep));
            tripletStateDot2.push_back(T(4 * (nStep - 1) + 4 * (i - 1) + j, 4 * i + j, -std::min(-control(i * 2 + 1) + yt * control(i * 2), 0.0)*100));
            tripletStateDot2.push_back(T(4 * 2 * (nStep - 1) + 4 * (i - 1) + j, 4*i+j, -std::min(control(i * 2 + 1) - yt * control(i * 2), 0.0)*100));
          }
          else if (var_set == "stateDot3")
          {
            tripletStateDot3.push_back(T(4 * (i - 1) + j, 4 * i + j, -tStep));
            tripletStateDot3.push_back(T(4 * 3 * (nStep - 1) + 4 * (i - 1) + j, 4*i+j, -std::min(control(i * 2 + 1) - yb * control(i * 2), 0.0)*100));
            tripletStateDot3.push_back(T(4 * 4 * (nStep - 1) + 4 * (i - 1) + j, 4*i+j, -std::min(-control(i * 2 + 1) + yb * control(i * 2), 0.0)*100));
          }
        }
        if (var_set == "state")
        {
          // df1/dtheta
          tripletState.push_back(T(4 * (i - 1) + 0, 4 * i + 2, - dStickdtheta(0)*tStep)); 
          tripletState.push_back(T(4 * (i - 1) + 1, 4 * i + 2, - dStickdtheta(1)*tStep)); 
          // df2/dtheta
          tripletState.push_back(T(4 * (nStep - 1) + 4 * (i - 1) + 0, 4 * i + 2, -control(i*2+1)+yt*control(i*2) < 0? (-control(i * 2 + 1) + yt * control(i * 2)) * dUpdtheta(0)*100 : 0.0));
          tripletState.push_back(T(4 * (nStep - 1) + 4 * (i - 1) + 1, 4 * i + 2, -control(i*2+1)+yt*control(i*2) < 0? (-control(i * 2 + 1) + yt * control(i * 2)) * dUpdtheta(1)*100 : 0.0));
          // df2/dtheta
          tripletState.push_back(T(T(4 * 3 * (nStep - 1) + 4 * (i - 1) + 0, 4*i+2, control(i * 2 + 1) - yb * control(i * 2) < 0? (control(i * 2 + 1) - yb * control(i * 2))*dDowndtheta(0)*100 : 0.0)));
          tripletState.push_back(T(T(4 * 3 * (nStep - 1) + 4 * (i - 1) + 1, 4*i+2, control(i * 2 + 1) - yb * control(i * 2) < 0? (control(i * 2 + 1) - yb * control(i * 2))*dDowndtheta(1)*100 : 0.0)));
        }
        if (var_set == "state") {
          jac_block.setFromTriplets(tripletState.begin(), tripletState.end());
        } else if (var_set == "control") 
        {
          jac_block.setFromTriplets(tripletControl.begin(), tripletControl.end());
        } else if (var_set == "stateDot2") 
        {
          jac_block.setFromTriplets(tripletStateDot2.begin(), tripletStateDot2.end());
        } else if (var_set == "stateDot3") 
        {
          jac_block.setFromTriplets(tripletStateDot3.begin(), tripletStateDot3.end());
        }
      }
    }

  private:
    double px, fMax, mMax, c, mu, m, muGround, tStep;
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
      YAML::Node params = YAML::LoadFile("/home/mzwang/qsp_ws/src/pusher/Config/params.yaml");
      QFinal = params["QFinal"].as<double>();
      Q = params["Q"].as<double>();
      R = params["R"].as<double>();
      weightQ.resize(4, 4);
      weightR.resize(2, 2);
      weightQ = Eigen::Vector4d(1, 3, 0.1, 0).asDiagonal();
      weightR = Eigen::Vector2d(1, 1).asDiagonal();
    }

    double GetCost() const override
    {
      VectorXd state = GetVariables()->GetComponent("state")->GetValues();
      VectorXd control = GetVariables()->GetComponent("control")->GetValues();
      VectorXd stateDot2 = GetVariables()->GetComponent("stateDot2")->GetValues();
      VectorXd stateDot3 = GetVariables()->GetComponent("stateDot3")->GetValues();
      int nStep = control.size() / 2;
      double cost = 0;
      for (int i = 0; i < nStep - 1; ++i)
      {
        cost += (state.segment(i * 4, 4) - stateNominal.segment(i * 4, 4)).transpose() *
                Q * weightQ * (state.segment(i * 4, 4) - stateNominal.segment(i * 4, 4));
        cost += (control.segment(i * 2, 2) - controlNominal.segment(i * 2, 2)).transpose() *
                R * weightR * (control.segment(i * 2, 2) - controlNominal.segment(i * 2, 2));
      }
      cost += (state.segment((nStep - 1) * 4, 4) - stateNominal.segment((nStep - 1) * 4, 4)).transpose() *
              QFinal * weightQ * (state.segment((nStep - 1) * 4, 4) - stateNominal.segment((nStep - 1) * 4, 4));

      // cost += 10.0*(stateDot2.lpNorm<1>() + stateDot3.lpNorm<1>());
      return cost;
    }

    void FillJacobianBlock(std::string var_set, Jacobian &jac) const override
    {
      VectorXd state = GetVariables()->GetComponent("state")->GetValues();
      VectorXd control = GetVariables()->GetComponent("control")->GetValues();
      std::vector<T> tripletState, tripletControl;
      int nStep = control.size() / 2;

      if (var_set == "state")
      {
        for (int i = 0; i < nStep - 1; ++i)
        {
          tripletState.push_back(T(0, i*4, 2*Q*weightQ(0,0)*(state(i*4)-stateNominal(i*4))));
          tripletState.push_back(T(0, i*4+1, 2*Q*weightQ(1,1)*(state(i*4+1)-stateNominal(i*4+1))));
          tripletState.push_back(T(0, i*4+2, 2*Q*weightQ(2,2)*(state(i*4+2)-stateNominal(i*4+2))));
          tripletState.push_back(T(0, i*4+3, 2*Q*weightQ(3,3)*(state(i*4+3)-stateNominal(i*4+3))));
        }
        tripletState.push_back(T(0, (nStep-1)*4, 2*QFinal*weightQ(0,0)*(state((nStep-1)*4)-stateNominal((nStep-1)*4))));
        tripletState.push_back(T(0, (nStep-1)*4+1, 2*QFinal*weightQ(1,1)*(state((nStep-1)*4+1)-stateNominal((nStep-1)*4+1))));
        tripletState.push_back(T(0, (nStep-1)*4+2, 2*QFinal*weightQ(2,2)*(state((nStep-1)*4+2)-stateNominal((nStep-1)*4+2))));
        tripletState.push_back(T(0, (nStep-1)*4+3, 2*QFinal*weightQ(3,3)*(state((nStep-1)*4+3)-stateNominal((nStep-1)*4+3))));
        jac.setFromTriplets(tripletState.begin(), tripletState.end());
      }
      else if (var_set == "control")
      {
        for (int i = 0; i < nStep - 1; ++i)
        {
          tripletControl.push_back(T(0, i*2, 2*R*weightR(0,0)*(control(i*2) - controlNominal(i*2))));
          tripletControl.push_back(T(0, i*2+1, 2*R*weightR(0,0)*(control(i*2+1) - controlNominal(i*2+1))));
        }
        jac.setFromTriplets(tripletControl.begin(), tripletControl.end());
      }
    }

  private:
    Eigen::VectorXd stateNominal, controlNominal;
    Eigen::MatrixXd weightQ, weightR;
    double QFinal, Q, R;
  };

} // namespace ifopt