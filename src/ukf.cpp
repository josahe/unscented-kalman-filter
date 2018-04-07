#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3; //TODO REVISIT

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5*M_PI; //TODO REVISIT

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  n_x_ = 5; //state dimension
  n_z_ = 3; //radar measurement dimension
  n_aug_ = 7; //augmented state dimension
  lambda_ = 3 - n_aug_; //spreading parameter
  n_sig_ = 2 * n_aug_ + 1; // number of sigma points

  //vector containing weights of sigma points
  weights_ = VectorXd(n_sig_);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < (n_sig_); ++i) {
    weights_(i) = 0.5 / (lambda_ + n_aug_);
  }

  //matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd(n_x_, n_sig_);

  //state covariance matrix TODO REVISIT using identity matrix
  P_.setIdentity();

  // access to normalise function for angles
  Tools tools;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  /*****************************************************************************
   Initialisation
  *****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    cout << "UKF: Intialised with ";

    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      cout << "LIDAR" << endl;
      //Initialize state
      x_ << meas_package.raw_measurements_[0], //px
            meas_package.raw_measurements_[1], //py
            0.0, //v
            0.0, //yaw
            0.0; //yawdot TODO REVISIT initialisation values
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      cout << "RADAR" << endl;
      //Convert radar from polar to cartesian coordinates and initialize state
      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
      x_ << rho * cos(phi), //px
            rho * sin(phi), //py
            0.0, //v
            0.0, //yaw
            0.0; //yawdot TODO REVISIT initialisation values
    }

    time_us_ = meas_package.timestamp_; //update previous timestamp
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   Prediction
  *****************************************************************************/
  //compute time elapsed between current and previous measurements
  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  Prediction(dt);

  /*****************************************************************************
   Update
  *****************************************************************************/
  if (use_laser_ && meas_package.sensor_type_==MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  }
  else if (use_radar_ && meas_package.sensor_type_==MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  //sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_);
  Xsig_aug.fill(0.0);

  //Generate augmented sigma points
  GenerateAugmentedSigmaPoints(&Xsig_aug);

  //Predict sigma points
  SigmaPointPrediction(Xsig_aug, delta_t);

  //Predict state mean and covariance
  PredictMeanAndCovariance();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement
 * and Kalman Filter equations
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  //new measurement
  VectorXd z = VectorXd(5);
  z << meas_package.raw_measurements_[0], //px
       meas_package.raw_measurements_[1], //py
       0.0, //v
       0.0, //yaw
       0.0; //yawdot TODO REVISIT initialisation values

  //measurement matrix
  MatrixXd H = MatrixXd(2, 5);

  //measurement covariance matrix
  MatrixXd R = MatrixXd(2, 2);

  // initializing matrices
  H << 1, 0, 0, 0, 0,
       0, 1, 0, 0, 0;
  R << 0.0225, 0,
       0, 0.0225;

  VectorXd z_pred = H * x_;
  VectorXd y = z.head(2) - z_pred;
  MatrixXd Ht = H.transpose();
  MatrixXd PHt = P_ * Ht; 
  MatrixXd S = H * PHt + R;
  MatrixXd Si = S.inverse();
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);

  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H) * P_;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

  //new measurement
  VectorXd z = VectorXd(n_z_);
  float rho = meas_package.raw_measurements_[0];
  float phi = meas_package.raw_measurements_[1];
  float rhodot = meas_package.raw_measurements_[2];
  z << rho,
       phi,
       rhodot;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_, n_sig_);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);

  //measurement covariance matrix S
  MatrixXd S_pred = MatrixXd(n_z_, n_z_);

  PredictRadarMeasurement(&Zsig, &z_pred, &S_pred);

  UpdateState(z, Zsig, z_pred, S_pred);
}

/**
 * Generates a matrix of sigma points, augmented with process noise modelling
 * @param {MatrixXd*} Xsig_aug
 */
void UKF::GenerateAugmentedSigmaPoints(MatrixXd* Xsig_aug) {

  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);
  x_aug.fill(0.0);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);
  P_aug.fill(0.0);

  //create augmented mean state
  x_aug.head(5) = x_;

  //create augmented covariance matrix
  MatrixXd Q = MatrixXd(2, 2); // process noise covariance matrix
  Q << std_a_*std_a_, 0.0,
       0.0, std_yawdd_*std_yawdd_;

  P_aug.topLeftCorner(5, 5) = P_;
  P_aug.bottomRightCorner(2, 2) = Q;

  //create square root matrix
  MatrixXd A = MatrixXd(7, 7);
  A = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug->col(0) = x_aug;

  for (int i = 0; i < n_aug_; ++i) {
    Xsig_aug->col(i + 1)          << x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
    Xsig_aug->col(i + 1 + n_aug_) << x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
  }
}

/**
 * Predicts a matrix of sigma points using the CTRV process model
 * @param {MatrixXd} Xsig_aug
 * @param {double} delta_t
 */
void UKF::SigmaPointPrediction(MatrixXd Xsig_aug, double delta_t) {

  Xsig_pred_.fill(0.0);

  //predict sigma points
  for (int i = 0; i < (n_sig_); ++i) {
    //extract state variables for better readability
    // TODO can define these variables as constants so compiler can optimise
    double px = Xsig_aug(0, i);
    double py = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawdot = Xsig_aug(4, i);
    double nu_acc = Xsig_aug(5, i);
    double nu_yaw = Xsig_aug(6, i);

    //predicted state variables
    double px_pred, py_pred;

    //avoid division by zero
    if (fabs(yawdot) > 0.001) {
      px_pred = px + v/yawdot * ( sin(yaw + yawdot*delta_t) - sin(yaw));
      py_pred = py + v/yawdot * (-cos(yaw + yawdot*delta_t) + cos(yaw));
    } else {
      px_pred = px + v * cos(yaw) * delta_t;
      py_pred = py + v * sin(yaw) * delta_t;
    }

    double v_pred = v;
    double yaw_pred = yaw + yawdot * delta_t;
    double yawdot_pred = yawdot;

    //add process noise
    px_pred += (0.5 * pow(delta_t, 2) * cos(yaw) * nu_acc);
    py_pred += (0.5 * pow(delta_t, 2) * sin(yaw) * nu_acc);
    v_pred += (delta_t * nu_acc);
    yaw_pred += (0.5 * pow(delta_t, 2) * nu_yaw);
    yawdot_pred += (delta_t * nu_yaw);

    //write predicted sigma points into right column
    Xsig_pred_.col(i) << px_pred,
                         py_pred,
                         v_pred,
                         yaw_pred,
                         yawdot_pred;
  }
}

/**
 * Predicts the state mean and covariance matrices
 */
void UKF::PredictMeanAndCovariance() {

  //predict state mean
  x_.fill(0.0);
  for (int i = 0; i < (n_sig_); ++i) {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  //predict state covariance matrix
  VectorXd X_x = VectorXd(n_x_);
  P_.fill(0.0);
  for (int i = 0; i < (n_sig_); ++i) {
    //state difference
    X_x = Xsig_pred_.col(i) - x_;

    //angle normalisation
    while (X_x(3) > M_PI) {X_x(3) -= 2.*M_PI;}
    while (X_x(3) <-M_PI) {X_x(3) += 2.*M_PI;}

    //state covariance
    P_ += weights_(i) * X_x * X_x.transpose();
  }
}

/**
 * Predicts radar measurements
 * @param {VectorXd*} Zsig
 * @param {MatrixXd*} z_pred
 * @param {MatrixXd*} S_pred
 */
void UKF::PredictRadarMeasurement(MatrixXd* Zsig, VectorXd* z_pred,
                                  MatrixXd* S_pred) {

  //measurement noise covariance matrix R
  MatrixXd R = MatrixXd(n_z_, n_z_);
  R.fill(0.0);
  R(0,0) = pow(std_radr_, 2);
  R(1,1) = pow(std_radphi_, 2);
  R(2,2) = pow(std_radrd_, 2);

  //transform sigma points into measurement space
  for (int i = 0; i < (n_sig_); ++i) {
    //extract state variables for better readability
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);
    //double yawdot = Xsig_pred(4, i);

    //predicted radar measurement variables
    double rho,    //radial distance
           phi,    //angle
           rhodot; //radial velocity

    rho = sqrt(pow(px, 2) + pow(py, 2)); // TODO atan2(0.0,0.0) is undefined
    phi = atan2(py, px);

    //angle normalisation
    while (phi > M_PI) {phi -= 2.*M_PI;}
    while (phi <-M_PI) {phi += 2.*M_PI;}

    if (rho > 0.001) {
      rhodot = (px*cos(yaw)*v + py*sin(yaw)*v) / rho;
    } else {
      cout << "WARNING: Divide by zero prevented" << endl;
      rhodot = 0.1; //TODO what to do here?
    }

    Zsig->col(i) << rho,
                    phi,
                    rhodot;
  }

  //calculate mean predicted measurement
  z_pred->fill(0.0);
  for (int i = 0; i < (n_sig_); ++i) {
    *z_pred += weights_(i) * Zsig->col(i);
  }

  //calculate innovation covariance matrix S
  VectorXd Z_z = VectorXd(n_z_);
  S_pred->fill(0.0);
  for (int i = 0; i < (n_sig_); ++i) {
    //state difference
    Z_z = Zsig->col(i) - *z_pred;

    //angle normalisation
    while (Z_z(1) > M_PI) {Z_z(1) -= 2.*M_PI;}
    while (Z_z(1) <-M_PI) {Z_z(1) += 2.*M_PI;}

    //state covariance
    *S_pred += weights_(i) * Z_z * Z_z.transpose();
  }

  //add measurement noise covariance matrix
  *S_pred += R;
}

/**
 * Updates state using sensor (radar/lidar) measurements
 * @param {VectorXd} z
 * @param {MatrixXd} Zsig
 * @param {VectorXd} z_pred
 * @param {MatrixXd} S_pred
 */
void UKF::UpdateState(VectorXd z, MatrixXd Zsig, VectorXd z_pred,
                      MatrixXd S_pred) {

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);

  //calculate cross correlation matrix
  VectorXd X_x = VectorXd(n_x_);
  VectorXd Z_z = VectorXd(n_x_);
  Tc.fill(0.0);
  for (int i = 0; i < (n_sig_); ++i) {
    //state difference
    X_x = Xsig_pred_.col(i) - x_;
    Z_z = Zsig.col(i) - z_pred;

    //angle normalisation
    while (X_x(3) > M_PI) {X_x(3) -= 2.*M_PI;}
    while (X_x(3) <-M_PI) {X_x(3) += 2.*M_PI;}
    while (Z_z(1) > M_PI) {Z_z(1) -= 2.*M_PI;}
    while (Z_z(1) <-M_PI) {Z_z(1) += 2.*M_PI;}

    //state covariance
    Tc += weights_(i) * X_x * Z_z.transpose();
  }

  //calculate Kalman gain K
  MatrixXd K = MatrixXd(n_x_, n_z_);
  K = Tc * S_pred.inverse();

  //angle normalisation
  Z_z = z - z_pred;
  while (Z_z(1) > M_PI) {Z_z(1) -= 2.*M_PI;}
  while (Z_z(1) <-M_PI) {Z_z(1) += 2.*M_PI;}

  //update state mean and covariance matrix
  x_ = x_ + K*Z_z;
  P_ = P_ - K*S_pred*K.transpose();
}
