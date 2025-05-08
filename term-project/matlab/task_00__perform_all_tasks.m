close all
clear
clc

warning off

%% Perform All Tasks

% Create IMU Model
task_01__create_IMU_model

% Single-Axis Motion Experiment
task_02__run_single_axis_rate_table_motion_experiment
task_03__evaluate_SAM_gyro_experiment
task_04__evaluate_SAM_accel_experiment
task_05__evaluate_SAM_accel_experiment_v2

% % Multi-Axis Motion Experiment
% task_06__run_multi_axis_motion_experiment
% task_07__evaluate_MAM_gyro_experiment
% task_08__evaluate_MAM_accel_experiment