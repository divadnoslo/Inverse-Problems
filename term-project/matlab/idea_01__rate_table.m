close all
clear
clc


%% Idea 01: Rate Table Batched Estimation

% Constants
g = 9.81;  


%% Create IMU Model

% Instantiate Empty Class
imu = ImuModel();

% Set Accelerometer Parameters
imu.b_a = [0.01; 0.02; 0.03];
imu.M_a = [...
    0.001,   0.0001,  0.0002; ...
    0.0003,  0.002,   0.0004; ...
    0.0005,  0.0006,  0.003]; 

% Set Accel Noise Parameters
accelSigma = 1e-2;
imu.AccelWhiteNoiseOneSigmaX = accelSigma; 
imu.AccelWhiteNoiseOneSigmaY = accelSigma; 
imu.AccelWhiteNoiseOneSigmaZ = accelSigma; 

% Set Gyroscope Parameters
imu.b_g = [0.01; 0.02; 0.03];
imu.M_g = [...
    0.001,   0.0001,  0.0002; ...
    0.0003,  0.002,   0.0004; ...
    0.0005,  0.0006,  0.003]; 

% Set Gyro Noise Parameters
gyroSigma = 1e-2;
imu.GyroWhiteNoiseOneSigmaX = gyroSigma; 
imu.GyroWhiteNoiseOneSigmaY = gyroSigma; 
imu.GyroWhiteNoiseOneSigmaZ = gyroSigma; 


%% Create Traditional Rate Table Motion

ambientTemperature = 21;
icm = ImuCalibrationManager(...
    "SamplingRate", 100, ...
    "Duration", 10);

dataset = icm.createCalibrationDataset(imu, ambientTemperature);

[Ga, da] = dataset2accelInverseProblem(dataset);
[Gg, dg] = dataset2accelInverseProblem(dataset);

% Rank
pa = rank(Ga)
pg = rank(Gg)

% Obviously, we can't continue

