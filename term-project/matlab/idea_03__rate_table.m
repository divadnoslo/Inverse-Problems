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
imu.b_a = [...
    6.4e-3 * g; ...
    -5.2e-3 * g; ...
    3.70e-3 * g];
imu.M_a = [...
    150e-6,   0.0001,  -0.0002; ...
    0.0003,  -175e-6,   -0.0004; ...
    0.0005,  -0.0006,  198e-6]; 

% Set Accel Noise Parameters
VRW = 0.07;
accelSigma = VRW / 60;
imu.AccelWhiteNoiseOneSigmaX = accelSigma; 
imu.AccelWhiteNoiseOneSigmaY = accelSigma; 
imu.AccelWhiteNoiseOneSigmaZ = accelSigma; 

% Set Gyroscope Parameters
imu.b_g = [100; 30; -250] * pi/180 / 60^2;
imu.M_g = [...
    450e-6,   -0.0001,  0.0002; ...
    -0.0003,  -300e-6,   0.0004; ...
    -0.0005,  0.0006,  175e-6]; 

% Set Gyro Noise Parameters
ARW = 0.15 * pi/180;
gyroSigma = ARW / 60;
imu.GyroWhiteNoiseOneSigmaX = gyroSigma; 
imu.GyroWhiteNoiseOneSigmaY = gyroSigma; 
imu.GyroWhiteNoiseOneSigmaZ = gyroSigma; 

% Get True Accelerometer Model
m_accel_true = [...
    imu.AccelFixedBiasX; ...
    imu.AccelScaleFactorErrorX; ...
    imu.AccelMisalignmentXY; ...
    imu.AccelMisalignmentXZ; ...
    imu.AccelFixedBiasY; ...
    imu.AccelMisalignmentYX; ...
    imu.AccelScaleFactorErrorY; ...
    imu.AccelMisalignmentYZ; ...
    imu.AccelFixedBiasZ; ...
    imu.AccelMisalignmentZX; ...
    imu.AccelMisalignmentZY; ...
    imu.AccelScaleFactorErrorZ];

% Get True Gyroscope Model
m_accel_true = [...
    imu.AccelFixedBiasX; ...
    imu.AccelScaleFactorErrorX; ...
    imu.AccelMisalignmentXY; ...
    imu.AccelMisalignmentXZ; ...
    imu.AccelFixedBiasY; ...
    imu.AccelMisalignmentYX; ...
    imu.AccelScaleFactorErrorY; ...
    imu.AccelMisalignmentYZ; ...
    imu.AccelFixedBiasZ; ...
    imu.AccelMisalignmentZX; ...
    imu.AccelMisalignmentZY; ...
    imu.AccelScaleFactorErrorZ];


%% Create Rate Table Motion

g = [0; 0; -9.81];

Fs = 100;
dt = 1 / Fs;
t_segment = 0 : dt : 1;
K = length(t_segment);

A = (2*pi) * 45 * pi/180;
curve = A * sin(2*pi*t_segment);

w_b__i_b_true = zeros(3, 3*K);
w_b__i_b_true(1,1:K) = curve;
w_b__i_b_true(2,K + (1:K)) = curve;
w_b__i_b_true(3,2*K + (1:K)) = curve;

t = (0:length(w_b__i_b_true)-1) * dt;

euler = cumtrapz(t, w_b__i_b_true, 2);
dcm = euler2dcm(euler(1,:), euler(2,:), euler(3,:));

f_b__i_b_true = squeeze(pagemtimes(dcm, g));

% Vizualize Angular Velocity
fig = figure("Name", "Angular Velocity Profile");
tl = tiledlayout(3, 1, "Parent", fig);
title(tl, "Angular Velocity Profile")
ax = nexttile(1);
hold(ax, "on")
plot(t, 180/pi * w_b__i_b_true(1,:), 'r')
title("\omega_x")
xlabel("Time [sec]")
ylabel("[deg/sec]")
xlim([t(1) t(end)])
grid on
grid minor
ax = nexttile(2);
hold(ax, "on")
plot(t, 180/pi * w_b__i_b_true(2,:), 'g')
title("\omega_y")
xlabel("Time [sec]")
ylabel("[deg/sec]")
xlim([t(1) t(end)])
grid on
grid minor
ax = nexttile(3);
hold(ax, "on")
plot(t, 180/pi * w_b__i_b_true(3,:), 'b')
title("\omega_z")
xlabel("Time [sec]")
ylabel("[deg/sec]")
xlim([t(1) t(end)])
grid on
grid minor
linkaxes(tl.Children, 'x')

% Vizualize Euler Angle Profile
fig = figure("Name", "Euler Angle Profile");
tl = tiledlayout(3, 1, "Parent", fig);
title(tl, "Euler Angle Profile")
ax = nexttile(1);
hold(ax, "on")
plot(t, 180/pi * euler(1,:), 'r')
title("Roll (\phi)")
xlabel("Time [sec]")
ylabel("[deg]")
xlim([t(1) t(end)])
grid on
grid minor
ax = nexttile(2);
hold(ax, "on")
plot(t, 180/pi * euler(2,:), 'g')
title("Pitch (\theta)")
xlabel("Time [sec]")
ylabel("[deg]")
xlim([t(1) t(end)])
grid on
grid minor
ax = nexttile(3);
hold(ax, "on")
plot(t, 180/pi * euler(3,:), 'b')
title("Yaw (\psi)")
xlabel("Time [sec]")
ylabel("[deg]")
xlim([t(1) t(end)])
grid on
grid minor
linkaxes(tl.Children, 'x')

% Vizualize Specific Force Profile
fig = figure("Name", "Specific Force Profile");
tl = tiledlayout(3, 1, "Parent", fig);
title(tl, "Specific Force Profile")
ax = nexttile(1);
hold(ax, "on")
plot(t, f_b__i_b_true(1,:), 'r')
title("f_x")
xlabel("Time [sec]")
ylabel("[m/sec/sec]")
xlim([t(1) t(end)])
grid on
grid minor
ax = nexttile(2);
hold(ax, "on")
plot(t, f_b__i_b_true(2,:), 'g')
title("f_y")
xlabel("Time [sec]")
ylabel("[m/sec/sec]")
xlim([t(1) t(end)])
grid on
grid minor
ax = nexttile(3);
hold(ax, "on")
plot(t, f_b__i_b_true(3,:), 'b')
title("f_z")
xlabel("Time [sec]")
ylabel("[m/sec/sec]")
xlim([t(1) t(end)])
grid on
grid minor
linkaxes(tl.Children, 'x')


%% Create IMU Measurements

[f_b__i_b_meas, w_b__i_b_meas] = imu.runForwardModel(f_b__i_b_true, w_b__i_b_true, 21);


%% Create Accelerometer Inverse Problem

% Create Model Operator
F = [ones(length(f_b__i_b_true), 1), f_b__i_b_true.'];
Z = zeros(size(F));
G = [F, Z, Z; Z, F, Z; Z, Z, F];
[m, n] = size(G);

% Create Data
d = [...
    f_b__i_b_meas(1,:).' - f_b__i_b_true(1,:).'; ...
    f_b__i_b_meas(2,:).' - f_b__i_b_true(2,:).'; ...
    f_b__i_b_meas(3,:).' - f_b__i_b_true(3,:).'];

% Make a Weighted Least Squares Problem
W = accelSigma * eye(m);
Gw = W * G;
dw = W * d;


%% L2 Regression

% L2 Regression
m_L2 = (Gw.' * Gw) \ (Gw.' * dw);

% L2 Model Regression Error
m_error = m_L2 - m_accel_true;

% Model Covariance
C = inv(Gw.' * Gw);

% 95% Confidence Bounds
conf95 = 1.96 * sqrt(diag(C));

% Print Results to Table
results = array2table([...
    m_accel_true, ...
    m_L2, ...
    conf95, ...
    m_error]);
results.Properties.VariableNames = {...
    'TrueModel', ...
    'L2RegressionModel', ...
    'Confidence95Bound', ...
    'L2ModelError'};
results.Properties.RowNames = {...
    'AccelFixedBiasX', ...
    'AccelScaleFactorErrorX', ...
    'AccelMisalignmentXY', ...
    'AccelMisalignmentXZ', ...
    'AccelFixedBiasY', ...
    'AccelMisalignmentYX', ...
    'AccelScaleFactorErrorY', ...
    'AccelMisalignmentYZ', ...
    'AccelFixedBiasZ', ...
    'AccelMisalignmentZX', ...
    'AccelMisalignmentZY', ...
    'AccelScaleFactorErrorZ'};
fprintf("Results\n\n")
disp(results)


%% SVD

% Is G full rank?
p = rank(G);
fprintf("Rank of G: %d\n\n", p)

% SVD
[U, S, V] = svd(G, 'econ'); 

% Truncate
Up = U(:,1:p);
Sp = S(1:p,1:p);
Vp = V(:,1:p);

% Generalized Inverse Solution
m_svd = Vp * inv(Sp) * Up.' * d;