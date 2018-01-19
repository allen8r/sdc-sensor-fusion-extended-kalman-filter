#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  if (estimations.size() == 0) {
    cout << "Error: estimation vector is empty!" << endl;
    return rmse;
  }

  if (estimations.size() != ground_truth.size()) {
    cout << "Error: estimation vector size does not equal ground truth vector size!" << endl;
    return rmse;
  }

  for (int i = 0; i < estimations.size(); ++i) {
    VectorXd residual = estimations[i] - ground_truth[i];
    // <vector>.array() allows element-wise operations
    VectorXd residual_squared = residual.array() * residual.array();
    rmse += residual_squared;
  }

  VectorXd mean_squared_error = rmse / estimations.size();
  rmse = mean_squared_error.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
}
