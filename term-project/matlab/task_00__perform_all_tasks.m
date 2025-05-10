close all
clear
clc

warning off

%% Perform All Tasks

% Create IMU Model
task_01__create_IMU_model

% Motion Profiles for Comparison
task_02__run_motion_profile_1
task_03__run_motion_profile_2
task_04__run_motion_profile_3
task_05__compare_gyroscope_results
task_06__compare_accelerometer_results

% Poorly-Scaled Case
task_07__run_motion_profile_4