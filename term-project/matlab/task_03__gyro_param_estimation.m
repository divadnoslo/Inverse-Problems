close all
clear
clc

% Save figures as *.eps
saveFigureAsEps = @(name, fig)(exportgraphics(fig, fullfile("..", "latex", "images", name)));

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
m_gyro_true = [...
    imu.GyroFixedBiasX; ...
    imu.GyroScaleFactorErrorX; ...
    imu.GyroMisalignmentXY; ...
    imu.GyroMisalignmentXZ; ...
    imu.GyroFixedBiasY; ...
    imu.GyroMisalignmentYX; ...
    imu.GyroScaleFactorErrorY; ...
    imu.GyroMisalignmentYZ; ...
    imu.GyroFixedBiasZ; ...
    imu.GyroMisalignmentZX; ...
    imu.GyroMisalignmentZY; ...
    imu.GyroScaleFactorErrorZ];


%% Create Rate Table Motion

g = [0; 0; -9.81];

Fs = 100;
dt = 1 / Fs;
t_segment = 0 : dt : 1;
K = length(t_segment);

A = (2*pi) * 15 * pi/180;
curve = A * sin(2*pi*t_segment);

w_b__i_b_true = zeros(3, 3*K);
w_b__i_b_true(1,1:K) = curve;
w_b__i_b_true(2,K + (1:K)) = curve;
w_b__i_b_true(3,2*K + (1:K)) = curve;

t = (0:length(w_b__i_b_true)-1) * dt;

euler = cumtrapz(t, w_b__i_b_true, 2);
dcm = euler2dcm(euler(1,:), euler(2,:), euler(3,:));

f_b__i_b_true = squeeze(pagemtimes(dcm, g));


%% Create IMU Measurements

[f_b__i_b_meas, w_b__i_b_meas] = imu.runForwardModel(f_b__i_b_true, w_b__i_b_true, 21);


%% Vizualize

% Vizualize Angular Velocity
fig = figure("Name", "Angular Velocity Profile");
tl = tiledlayout(3, 1, "Parent", fig);
title(tl, "Angular Velocity Profile")
ax = nexttile(1);
hold(ax, "on")
plot(t, 180/pi * w_b__i_b_true(1,:), 'k', 'LineWidth', 2)
plot(t, 180/pi * w_b__i_b_meas(1,:), 'r')
title("\omega_x")
xlabel("Time [sec]")
ylabel("[deg/sec]")
xlim([t(1) t(end)])
grid on
grid minor
legend(["Test Bed", "UUT"], "Location", "eastoutside")
ax = nexttile(2);
hold(ax, "on")
plot(t, 180/pi * w_b__i_b_true(2,:), 'k', 'LineWidth', 2)
plot(t, 180/pi * w_b__i_b_meas(2,:), 'g')
title("\omega_y")
xlabel("Time [sec]")
ylabel("[deg/sec]")
xlim([t(1) t(end)])
grid on
grid minor
legend(["Test Bed", "UUT"], "Location", "eastoutside")
ax = nexttile(3);
hold(ax, "on")
plot(t, 180/pi * w_b__i_b_true(3,:), 'k', 'LineWidth', 2)
plot(t, 180/pi * w_b__i_b_meas(3,:), 'b')
title("\omega_z")
xlabel("Time [sec]")
ylabel("[deg/sec]")
xlim([t(1) t(end)])
grid on
grid minor
legend(["Test Bed", "UUT"], "Location", "eastoutside")
linkaxes(tl.Children, 'x')
saveFigureAsEps("SAM_angular_velocity_profile.eps", fig)

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
saveFigureAsEps("SAM_euler_angle_profile.eps", fig)

% Vizualize Specific Force Profile
fig = figure("Name", "Specific Force Profile");
tl = tiledlayout(3, 1, "Parent", fig);
title(tl, "Specific Force Profile")
ax = nexttile(1);
hold(ax, "on")
plot(t, f_b__i_b_true(1,:), 'k', 'LineWidth', 2)
plot(t, f_b__i_b_meas(1,:), 'r')
title("f_x")
xlabel("Time [sec]")
ylabel("[m/sec/sec]")
xlim([t(1) t(end)])
grid on
grid minor
legend(["Test Bed", "UUT"], "Location", "eastoutside")
ax = nexttile(2);
hold(ax, "on")
plot(t, f_b__i_b_true(2,:), 'k', 'LineWidth', 2)
plot(t, f_b__i_b_meas(2,:), 'g')
title("f_y")
xlabel("Time [sec]")
ylabel("[m/sec/sec]")
xlim([t(1) t(end)])
grid on
grid minor
legend(["Test Bed", "UUT"], "Location", "eastoutside")
ax = nexttile(3);
hold(ax, "on")
plot(t, f_b__i_b_true(3,:), 'k', 'LineWidth', 2)
plot(t, f_b__i_b_meas(3,:), 'b')
title("f_z")
xlabel("Time [sec]")
ylabel("[m/sec/sec]")
xlim([t(1) t(end)])
grid on
grid minor
legend(["Test Bed", "UUT"], "Location", "eastoutside")
linkaxes(tl.Children, 'x')
saveFigureAsEps("SAM_specific_force_profile.eps", fig)


%% Create Gyroscope Inverse Problem

% Create Model Operator
Omega = [ones(length(w_b__i_b_true), 1), w_b__i_b_true.'];
Z = zeros(size(Omega));
G = [Omega, Z, Z; Z, Omega, Z; Z, Z, Omega];
[m, n] = size(G);

% Create Data
d = [...
    w_b__i_b_meas(1,:).' - w_b__i_b_true(1,:).'; ...
    w_b__i_b_meas(2,:).' - w_b__i_b_true(2,:).'; ...
    w_b__i_b_meas(3,:).' - w_b__i_b_true(3,:).'];

% Make a Weighted Least Squares Problem
W = gyroSigma * eye(m);
Gw = W * G;
dw = W * d;


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
m_error = m_svd - m_gyro_true;

% Model Resolution
Rm = Vp.' * Vp;
diag(Rm);

% Data Resolution
Rd = Up.' * Up;
diag(Rd);

% Print Results to Table
results = array2table([...
    m_gyro_true, ...
    m_svd, ...
    m_error]);
results.Properties.VariableNames = {...
    'TrueModel', ...
    'ModelParameters', ...
    'ModelError'};
results.Properties.RowNames = {...
    'GyroFixedBiasX', ...
    'GyroScaleFactorErrorX', ...
    'GyroMisalignmentXY', ...
    'GyroMisalignmentXZ', ...
    'GyroFixedBiasY', ...
    'GyroMisalignmentYX', ...
    'GyroScaleFactorErrorY', ...
    'GyroMisalignmentYZ', ...
    'GyroFixedBiasZ', ...
    'GyroMisalignmentZX', ...
    'GyroMisalignmentZY', ...
    'GyroScaleFactorErrorZ'};

% Plot Singular Values
fig = figure("Name", "SVD Singular Values");
ax = gca;
hold(ax, "on")
semilogy(1:n, diag(S), 'bo')
title("Gyroscope Singular Values")
xlabel("Singular Values")
ylabel("s_i")
grid on
grid minor
saveFigureAsEps("SAM_gyro_singular_values.eps", fig)


%% Display Results

fprintf("Results\n\n")
disp(results)

% Bar Chart Prep
params = categorical(results.Properties.RowNames);
params = reordercats(params, results.Properties.RowNames);
b = [1, 5, 9];
sf = [2, 7, 12];
ma = [3, 4, 6, 8, 10, 11];

% Bar Chart
fig = figure("Name", "Model Parameter Comparison");
tl = tiledlayout(3, 1, "Parent", fig);
title(tl, "Gyroscope Calibration Model Parameter Comparison")

ax = nexttile(1);
hold(ax, "on")
bh = bar(...
    ax, ...
    removecats(params(b)), ...
    180/pi * [results.TrueModel(b).'; results.ModelParameters(b).'].');
title("Fixed Biases")
ylabel("[deg/sec]")
grid on
grid minor
legend(["Truth", "Estimate"], "Location", "eastoutside")

ax = nexttile(2);
hold(ax, "on")
bar(...
    ax, ...
    removecats(params(sf)), ...
    1e6 * [results.TrueModel(sf).'; results.ModelParameters(sf).'].')
title("Scale Factor Errors")
ylabel("[ppm]")
grid on
grid minor
legend(["Truth", "Estimate"], "Location", "eastoutside")

ax = nexttile(3);
hold(ax, "on")
bar(...
    ax, ...
    removecats(params(ma)), ...
    1e3 * [results.TrueModel(ma).'; results.ModelParameters(ma).'].')
title("Misalignment Terms")
ylabel("[m-rad]")
grid on
grid minor
legend(["Truth", "Estimate"], "Location", "eastoutside")

saveFigureAsEps("SAM_gyro_model_parameters.eps", fig)

% Model Parameter Error
fig = figure("Name", "Model Parameter Error");
tl = tiledlayout(3, 1, "Parent", fig);
title(tl, "Gyroscope Parameter Error")

ax = nexttile(1);
hold(ax, "on")
bh = bar(...
    ax, ...
    removecats(params(b)), ...
    180/pi * results.ModelError(b));
title("Fixed Biases")
ylabel("[deg/sec]")
grid on
grid minor

ax = nexttile(2);
hold(ax, "on")
bar(...
    ax, ...
    removecats(params(sf)), ...
    1e6 * results.ModelError(sf))
title("Scale Factor Errors")
ylabel("[ppm]")
grid on
grid minor

ax = nexttile(3);
hold(ax, "on")
bar(...
    ax, ...
    removecats(params(ma)), ...
    1e3 * results.ModelError(ma))
title("Misalignment Terms")
ylabel("[m-rad]")
grid on
grid minor
saveFigureAsEps("SAM_gyro_model_error.eps", fig)


%% L2 Regression

% L2 Regression
m_L2 = (Gw.' * Gw) \ (Gw.' * dw);

% L2 Model Regression Error
m_error = m_L2 - m_gyro_true;

% Model Covariance
C = inv(Gw.' * Gw);

% 95% Confidence Bounds
conf95 = 1.96 * sqrt(diag(C));

% Vizualize
fig = figure("Name", "95 Confidence Bounds");
tl = tiledlayout(3, 1, "Parent", fig);
title(tl, "95%% Confidence Bounds")

ax = nexttile(1);
hold(ax, "on")
bh = bar(...
    ax, ...
    removecats(params(b)), ...
    180/pi * conf95(b));
title("Fixed Biases")
ylabel("[deg/sec]")
grid on
grid minor

ax = nexttile(2);
hold(ax, "on")
bar(...
    ax, ...
    removecats(params(sf)), ...
    1e6 * conf95(sf))
title("Scale Factor Errors")
ylabel("[ppm]")
grid on
grid minor

ax = nexttile(3);
hold(ax, "on")
bar(...
    ax, ...
    removecats(params(ma)), ...
    conf95(ma))
title("Misalignment Terms")
ylabel("[rad]")
grid on
grid minor
saveFigureAsEps("SAM_gyro_model_95_confidence_bounds.eps", fig)