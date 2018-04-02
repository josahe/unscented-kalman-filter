#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

double Tools::NormaliseAngle(double angle) {
  while (angle > M_PI) {angle -= 2.*M_PI;}
  while (angle <-M_PI) {angle += 2.*M_PI;}
  return angle;
}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() == 0) {
    cout << "Estimation vector is empty" << endl;
  } else if (estimations.size() != ground_truth.size()) {
    cout << "Vecor size does not match" << endl;
  }

  //accumulate squared residuals
  for (unsigned int i = 0; i < estimations.size(); ++i)  {
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse += residual;
  }

  //calculate the mean
  rmse = rmse.array() / estimations.size();

  //calculate the squared root
  rmse = rmse.array().sqrt();

  assert (~rmse.isZero());
  return rmse;
}
