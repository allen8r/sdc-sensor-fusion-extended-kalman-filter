#include "kalman_filter.h"
#include <math.h>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  cout << "KalmanFilter::Predict()..." << endl;
  // predict the state
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * Q_ * Ft + Q_;
  cout << "P_:\n" << P_ << endl;
}

void KalmanFilter::Update(const VectorXd &z) {
  //update the state by using Kalman Filter equations
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred; // error
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  // new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  // update the state by using Extended Kalman Filter equations
  
  // Pre-conditions:
  // H_ is already set to Jacobian matrix Hj
  // h_x is already set to h(x') function
  cout << "KalmanFilter::UpdateEKF()..." << endl;
  cout << "z.size(): " << z.size() << endl << endl;
  cout << "h_x.size(): " << h_x.size() << endl << endl;
  VectorXd y = z - h_x;
  cout << "y:\n" << y << endl << endl;
  
  // Normalize angle phi
  float phi = y[1];
  // while phi not in expected range (-PI, PI)
  // keep adding 
  while (!(phi > -M_PI && phi < M_PI)) {
    phi += M_2_PI;
  }
  y[1] = phi;

  cout << "After normalizing phi..." << endl;
  cout << "y:\n" << y << endl << endl;

  MatrixXd Ht = H_.transpose();
  cout << "Ht:\n" << Ht << endl;
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  // new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
