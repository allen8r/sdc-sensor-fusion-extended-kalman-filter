#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
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

  // Add up the squared residuals
  for (int i = 0; i < estimations.size(); ++i) {
    VectorXd residual = estimations[i] - ground_truth[i];
    // <vector>.array() allows element-wise operations
    VectorXd residual_squared = residual.array() * residual.array();
    rmse += residual_squared;
  }

  // Find the mean of the squared residuals
  VectorXd mean_squared_error = rmse / estimations.size();

  // Take the square root of the mean squared error
  rmse = mean_squared_error.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj(3, 4);

  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  float px2 = px * px;
  float py2 = py * py;
  float px2_py2 = px2 + py2;
  float rt_px2py2 = sqrt(px2_py2);
  float px2py2_to_3_half = px2_py2 * rt_px2py2;

  // Check division by zero
  if (fabs(px2_py2) < 0.0001) {
    cout << "Error: Division by zero!" << endl;
    return Hj;
  }

  // Compute Jacobian matrix
  Hj << px / rt_px2py2,  py / rt_px2py2, 0, 0,
         -py / px2_py2,    px / px2_py2, 0, 0,
          py * (vx*py - vy*px) / px2py2_to_3_half, px * (vy*px - vx*py) / px2py2_to_3_half, px / rt_px2py2, py / rt_px2py2;

  return Hj; 
}

VectorXd Tools::Calculate_h_x(const VectorXd& x_state) {
  VectorXd h_x = VectorXd(3);
  h_x << 0, 0, 0;
  
  float px = x_state[0];
  float py = x_state[1];
  float vx = x_state[2];
  float vy = x_state[3];

  float px2 = px * px;
  float py2 = py * py;
  float px2_py2 = px2 + py2;
  
  if (fabs(px) < 0.0001 || fabs(px2_py2) < 0.0001) {
    cout << "Error: division by zero in h(x)" << endl;
    return h_x;
  } 

  // Compute the function h(x')
  h_x << sqrt(px2_py2), atan2(py, px), (px*vx + py*vy) / sqrt(px2_py2);

  return h_x;
}
