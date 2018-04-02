# Unscented Kalman Filter
> A sensor data processing pipeline to estimate the state of a moving object

## INTRODUCTION
The goal of this project is to write an Unscented Kalman Filter to estimate the state of a moving object using noisy lidar and radar measurements. To predict the location of the object, a [Constant Turn Rate and Velocity (CTRV) model](https://www.researchgate.net/publication/4370048_Comparison_and_evaluation_of_advanced_motion_models_for_vehicle_tracking) is used.

This project was undertaken as part of the [Udacity Self-Driving Car NanoDegree](https://eu.udacity.com/course/self-driving-car-engineer-nanodegree--nd013).

### Pipeline summary
The Unscented Kalman Filter (UKF), like the [Extended Kalman Filter](https://github.com/josahe/extended-kalman-filter), is a continuous loop of measurement updates and predictions. The UKF deals with non-linearity differently to the EKF, by mapping selected points in the probability distribution graph (sigma points) through the non-linear measurement or non-linear process function, and then calculating the mean state vector and covariance matrix to produce a best-fit linear (Gaussian) representation of the location of the object.

The CTRV model predicts new (x,y) coordinates and turning angle using the assumption of a constant velocity and a constant turn rate (angular velocity).

## HOW TO USE
### Project dependencies
You can follow the guide in the README of the original project repo.
* https://github.com/udacity/CarND-Unscented-Kalman-Filter-Project

## RELEVANT FILES
* [ukf.cpp](src/ukf.cpp)
  * Initialises the state parameters and implements the UKF equations for prediction and radar measurements
