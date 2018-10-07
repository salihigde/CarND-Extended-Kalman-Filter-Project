#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
			  0, 1, 0, 0,
			  0, 0, 1000, 0,
			  0, 0, 0, 1000;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    previous_timestamp_ = measurement_pack.timestamp_;

    ekf_.Q_ = MatrixXd(4, 4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
     double x = measurement_pack.raw_measurements_[0] * cos(measurement_pack.raw_measurements_[1]);
     double y = measurement_pack.raw_measurements_[0] * sin(measurement_pack.raw_measurements_[1]);
     double vx = measurement_pack.raw_measurements_[2] * cos(measurement_pack.raw_measurements_[1]);
  	 double vy = measurement_pack.raw_measurements_[2] * sin(measurement_pack.raw_measurements_[1]);

     if(x <= 0){
       x = 0.00001;
     }
     if(y <= 0){
       y = 0.00001;
     }
     
     ekf_.x_ << x,y,vx,vy;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;

  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, dt, 0,
             0, 1, 0, dt,
             0, 0, 1, 0,
             0, 0, 0, 1;

	previous_timestamp_ = measurement_pack.timestamp_;

	float dt_2 = dt * dt;
	float dt_3 = dt_2 * dt;
	float dt_4 = dt_3 * dt;

  float noise_ax = 9.0;
  float noise_ay = 9.0;
  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
			   0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
			   dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
			   0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.H_ = Hj_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    // Radar updates
  } else {
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
