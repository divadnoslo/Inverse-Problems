close all
clear
clc

rng(0);


%% Demo 01: Demonstrate Basic IMU Model Use


%% Create IMU Model

% Constants
g = 9.81;  

% Instantiate Empty Class
imu = ImuModel();

% % Set Accelerometer Parameters
% sigma_a_fb = 10 * g * 1e-3;
% imu.b_a = sigma_a_fb * randn(3,1); 
% sigma_a_sf = 1000 * 1e-6;
% sigma_a_m = 0.05 * pi/180;
% imu.M_a = randn(3, 3) .* ...
%     [sigma_a_sf, sigma_a_m,  sigma_a_m; ...
%      sigma_a_m,  sigma_a_sf, sigma_a_m; ...
%      sigma_a_m,  sigma_a_m,  sigma_a_sf]; 
% sigma_a_wn = 0.0;
% imu.sigma_a = sigma_a_wn * ones(3, 1);
% 
% % Set Gyroscope Parameters
% sigma_g_fb = 1000 * pi/180 * (1/3600);
% imu.b_g = sigma_g_fb * randn(3,1); 
% sigma_g_sf = 1000 * 1e-6;
% sigma_g_m = 0.05 * pi/180;
% imu.M_g = randn(3, 3) .* ...
%     [sigma_g_sf, sigma_g_m,  sigma_g_m; ...
%      sigma_g_m,  sigma_g_sf, sigma_g_m; ...
%      sigma_g_m,  sigma_g_m,  sigma_g_sf]; 
% sigma_g_wn = 0 * pi/180;
% imu.sigma_g = sigma_g_wn * ones(3, 1);

% Set Accelerometer Parameters
imu.b_a = [0.01; 0.02; 0.03];
imu.M_a = [...
    0.001,   0.0001,  0.0002; ...
    0.0003,  0.002,   0.0004; ...
    0.0005,  0.0006,  0.003]; 

% Set Gyroscope Parameters
imu.b_g = [0.01; 0.02; 0.03];
imu.M_g = [...
    0.001,   0.0001,  0.0002; ...
    0.0003,  0.002,   0.0004; ...
    0.0005,  0.0006,  0.003];  


%% Create and Process Example Calibration Data

imuCalManager = ImuCalibrationManager();

calData = imuCalManager.createCalibrationDataSet(imu);
% imuCalManager.plotCalibrationDataset(calData);

processedImuModel = imuCalManager.processCalibrationDataSet(calData)