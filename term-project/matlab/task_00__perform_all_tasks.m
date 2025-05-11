close all
clear
clc

warning off

% Add "Lib" to Search Path
addpath(genpath(fullfile("..", "..", "PEIP-master", "Lib")))


%% Perform All Tasks

% Create IMU Model
task_01__create_IMU_model

% Traditional Calibration Comparison
task_02__run_traditional_motion_profile

% Motion Profiles for Comparison
task_03__run_motion_profile_1
task_04__run_motion_profile_2
task_05__run_motion_profile_3
task_06__compare_gyroscope_results
task_07__compare_accelerometer_results

% Poorly-Scaled Case
task_08__run_motion_profile_4