%% Task 08: Motion Profile 4

% close all
% clear
% clc

% Load Working File
load(fullfile(pwd, "working_file.mat"))

% Profile Number
profileNumber = 4;

% Helper Functions
saveFigureAsEps = @(name, fig)(exportgraphics(fig, fullfile("..", "latex", "images", name)));
preamble = @(profileNumber, description)(sprintf("Motion Profile %s: %s", num2str(profileNumber), description));
makeFileName = @(profileNumber, description)(sprintf("MP%s_%s", num2str(profileNumber), description));


%% Create Multi-Axis Rate Table Motion

g = [0; 0; -gravity];

Fs = 100;
dt = 1 / Fs;
t_end = 10;
t = 0 : dt : t_end;

curveX = ((2*pi) * 45 * pi/180) * sin(2*pi*(1/t_end)*t);
curveY = ((2*pi) * 30 * pi/180) * cos(2*pi*(1/t_end)*t);
curveZ = ((2*pi) * 15 * pi/180) * -sin(2*pi*(1/t_end)*t);

w_b__i_b_true = [...
    curveX, -curveX; ...
    curveY, -curveY; ...
    curveZ, -curveZ];
K = length(w_b__i_b_true);
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
title(tl, preamble(profileNumber, "Angular Velocity Profile"))
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
saveFigureAsEps(makeFileName(profileNumber, "angular_velocity_profile.eps"), fig)

% Vizualize Angular Velocity Error
fig = figure("Name", "Angular Velocity Error");
tl = tiledlayout(3, 1, "Parent", fig);
title(tl, preamble(profileNumber, "Angular Velocity Error"))
ax = nexttile(1);
hold(ax, "on")
plot(t, 180/pi * (w_b__i_b_meas(1,:) - w_b__i_b_true(1,:)), 'r', 'LineWidth', 2)
title("\Delta\omega_x")
xlabel("Time [sec]")
ylabel("[deg/sec]")
xlim([t(1) t(end)])
grid on
grid minor
ax = nexttile(2);
hold(ax, "on")
plot(t, 180/pi * (w_b__i_b_meas(2,:) - w_b__i_b_true(2,:)), 'r', 'LineWidth', 2)
title("\Delta\omega_y")
xlabel("Time [sec]")
ylabel("[deg/sec]")
xlim([t(1) t(end)])
grid on
grid minor
ax = nexttile(3);
hold(ax, "on")
plot(t, 180/pi * (w_b__i_b_meas(3,:) - w_b__i_b_true(3,:)), 'r', 'LineWidth', 2)
title("\Delta\omega_z")
xlabel("Time [sec]")
ylabel("[deg/sec]")
xlim([t(1) t(end)])
grid on
grid minor
linkaxes(tl.Children, 'x')
saveFigureAsEps(makeFileName(profileNumber, "angular_velocity_error.eps"), fig)

% Vizualize Euler Angle Profile
fig = figure("Name", "Euler Angle Profile");
tl = tiledlayout(3, 1, "Parent", fig);
title(tl, preamble(profileNumber, "Euler Angle Profile"))
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
saveFigureAsEps(makeFileName(profileNumber, "euler_angle_profile.eps"), fig)

% Vizualize Specific Force Profile
fig = figure("Name", "Specific Force Profile");
tl = tiledlayout(3, 1, "Parent", fig);
title(tl, preamble(profileNumber, "Specific Force Profile"))
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
saveFigureAsEps(makeFileName(profileNumber, "specific_force_profile.eps"), fig)

% Vizualize Specific Force Error
fig = figure("Name", "Specific Force Error");
tl = tiledlayout(3, 1, "Parent", fig);
title(tl, preamble(profileNumber, "Specific Force Error"))
ax = nexttile(1);
hold(ax, "on")
plot(t, f_b__i_b_meas(1,:) - f_b__i_b_true(1,:), 'r', 'LineWidth', 2)
title("\Delta f_x")
xlabel("Time [sec]")
ylabel("[m/sec/sec]")
xlim([t(1) t(end)])
grid on
grid minor
ax = nexttile(2);
hold(ax, "on")
plot(t, f_b__i_b_meas(2,:) - f_b__i_b_true(2,:), 'r', 'LineWidth', 2)
title("\Delta f_y")
xlabel("Time [sec]")
ylabel("[m/sec/sec]")
xlim([t(1) t(end)])
grid on
grid minor
ax = nexttile(3);
hold(ax, "on")
plot(t, f_b__i_b_meas(3,:) - f_b__i_b_true(3,:), 'r', 'LineWidth', 2)
title("\Delta f_z")
xlabel("Time [sec]")
ylabel("[m/sec/sec]")
xlim([t(1) t(end)])
grid on
grid minor
linkaxes(tl.Children, 'x')
saveFigureAsEps(makeFileName(profileNumber, "specific_force_error.eps"), fig)


%% Evaluate Gyro

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

% % Make This A Weighted Least Squares Problem
% W = accelSigma * eye(m);
% G = W * G;
% d = W * d;

% Examine Singular Values
[U, S, V] = svd(G, 'econ'); 
fig = figure("Name", "SVD Singular Values");
ax = gca;
hold(ax, "on")
ax.YScale = "log";
plot(1:n, diag(S), 'bo')
title(preamble(profileNumber, "Gyroscope Singular Values"))
xlabel("Singular Value Index")
ylabel("s_i")
grid on
grid minor
saveFigureAsEps(makeFileName(profileNumber, "gyro_singular_values.eps"), fig)

%% Try out Truncated SVD

% Truncated-SVD
p = 9;
Up = U(:, 1:p);
Sp = S(1:p, 1:p);
Vp = V(:, 1:p);
U0 = U(:, p+1:end);
S0 = S(p+1:end, p+1:end);
V0 = V(:, p+1:end);

% Compute SVD Solution
m_svd = Vp * inv(Sp) * Up.' * d;
gyroTable.SvdModel = m_svd;

% Model Error
gyroTable.SvdModelError = m_svd - m_gyro_true;

% Model Resolution
Rm = Vp * Vp.';

% Plot Model Resolution Matrix
fig = figure("Name", "Gyro Model Resolution Matrix");
ax = gca;
colormap('gray')
imagesc(Rm)
% clim([-1 1])
ax.XTick = 1 : n;
ax.YTick = 1 : n;
ax.XTickLabel = gyroModelLabels;
ax.YTickLabel = gyroModelLabels;
ax.YDir = "reverse";
axis equal
ax.XLim = [ax.XTick(1) - 0.5, ax.XTick(end) + 0.5];
ax.YLim = [ax.YTick(1) - 0.5, ax.YTick(end) + 0.5];
title(preamble(profileNumber, "Gyro Model Resolution Mtrix"))
colorbar(ax, "eastoutside")
saveFigureAsEps(makeFileName(profileNumber, "gyro_model_resolution.eps"), fig)

% Model Null Space
disp("Model Null Space")
disp(V0)


%% Tikohonov Regularization

s = diag(S);

% L-Curve
[rho, eta, reg_param] = l_curve_tikh_svd(U, s, d, 1000);
[alpha_tikh, ireg_corner] = l_curve_corner(rho, eta, reg_param);
rho_corner = rho(ireg_corner);
eta_corner = eta(ireg_corner);

% Report alpha_tikh
fprintf("alpha_tikh: %4.3e\n\n", alpha_tikh)

% Plot L-Curve
fig = figure("Name", "Zeroth-Order Tikhonov L-Curve");
ax = gca;
hold(ax, "on")
plot(ax, rho, eta, 'b', 'LineWidth', 3)
plot(ax, rho_corner, eta_corner, 'rx', 'MarkerSize', 20)
ax.XAxis.Scale = "log";
ax.YAxis.Scale = "log";
title(preamble(profileNumber, "Zeroth-Order Tikhonov L-Curve"))
xlabel("||Gm - d||_2")
ylabel("||m||_2")
grid on
grid minor
saveFigureAsEps(makeFileName(profileNumber, "accel_0th_order_tik_l_curve.eps"), fig)

% Compute Zeroth-Order Tikhonov Solution
m_tikh = (G.'*G + (alpha_tikh^2)*eye(n)) \ G.' * d;
gyroTable.TikhModel = m_tikh;

% Model Error
gyroTable.TikhModelError = m_tikh - m_gyro_true;


%% Examine Results

% Put Into One Struct
gyroResults(1) = mp1Gyro;
gyroResults(2) = mp2Gyro;
gyroResults(3) = mp3Gyro;

% Indexes
bIndexes = [1, 5, 9];
sfIndexes = [2, 7, 12];
mIndexes = [3, 4, 6, 8, 10, 11];

% Create Catagories
gyroBiases = categorical(gyroModelLabels(bIndexes));
gyroScaleFactors = categorical(gyroModelLabels(sfIndexes));
gyroMisalignments = categorical(gyroModelLabels(mIndexes));

% Set Colors
trueColor = "#000000";
svdColor = "#00FF00";
tikhColor = "#0000FF";

% Model Parameters
fig = figure("Name", "Model Parameter Comparison");
tl = tiledlayout(3, 1, "Parent", fig);
title(tl, preamble(profileNumber, "Gyro Parameter Estimation Comparison"))
ax = nexttile(1);
b = bar(...
    gyroBiases, ...
    [gyroTable.TrueModel(bIndexes).'; ...
     gyroTable.SvdModel(bIndexes).'; ...
     gyroTable.TikhModel(bIndexes).'].' * 180/pi);
b(1).FaceColor = trueColor;
b(2).FaceColor = svdColor;
b(3).FaceColor = tikhColor;
title("Gyroscope Bias Comparison")
ylabel("[deg/sec]")
grid on
grid minor
legend(["Truth", "SVD", "Tikh"], "Location", "eastoutside")
ax = nexttile(2);
b = bar(...
    gyroScaleFactors, ...
    [gyroTable.TrueModel(sfIndexes).'; ...
     gyroTable.SvdModel(sfIndexes).'; ...
     gyroTable.TikhModel(sfIndexes).'].' * 1e6);
b(1).FaceColor = trueColor;
b(2).FaceColor = svdColor;
b(3).FaceColor = tikhColor;
title("Gyroscope Scale Factor Error Comparison")
ylabel("[ppm]")
grid on
grid minor
legend(["Truth", "SVD", "Tikh"], "Location", "eastoutside")
ax = nexttile(3);
b = bar(...
    gyroMisalignments, ...
    [gyroTable.TrueModel(mIndexes).'; ...
     gyroTable.SvdModel(mIndexes).'; ...
     gyroTable.TikhModel(mIndexes).'].' * 1e3);
b(1).FaceColor = trueColor;
b(2).FaceColor = svdColor;
b(3).FaceColor = tikhColor;
title("Gyroscope Misalignment Comparison")
ylabel("[milli-rad]")
grid on
grid minor
legend(["Truth", "SVD", "Tikh"], "Location", "eastoutside")
saveFigureAsEps(makeFileName(profileNumber, "model_parameter_comparison.eps"), fig)

% Model Errors
fig = figure("Name", "Model Absolute Error");
tl = tiledlayout(3, 1, "Parent", fig);
title(tl, preamble(profileNumber, "Gyro Parameter Absolute Error"))
ax = nexttile(1);
b = bar(...
    gyroBiases, ...
    abs([gyroTable.SvdModelError(bIndexes).'; ...
         gyroTable.TikhModelError(bIndexes).']).' * 180/pi);
b(1).FaceColor = svdColor;
b(2).FaceColor = tikhColor;
title("Gyroscope Bias Absolute Error")
ylabel("[deg/sec]")
grid on
grid minor
legend(["SVD", "Tikh"], "Location", "eastoutside")
ax = nexttile(2);
b = bar(...
    gyroScaleFactors, ...
    abs([gyroTable.SvdModelError(sfIndexes).'; ...
         gyroTable.TikhModelError(sfIndexes).']) * 1e6);
ax.YScale = "log";
ax.YLim(1) = 1e-2;
b(1).FaceColor = svdColor;
b(2).FaceColor = tikhColor;
title("Gyroscope Scale Factor Absolute Error")
ylabel("[ppm]")
grid on
grid minor
legend(["SVD", "Tikh"], "Location", "eastoutside")
ax = nexttile(3);
b = bar(...
    gyroMisalignments, ...
    abs([gyroTable.SvdModelError(mIndexes).'; ...
     gyroTable.TikhModelError(mIndexes).']) * 1e3);
b(1).FaceColor = svdColor;
b(2).FaceColor = tikhColor;
title("Gyroscope Misalignment Absolute Error")
ylabel("[milli-rad]")
ax.YScale = "log";
ax.YLim(1) = 1e-7;
grid on
grid minor
legend(["SVD", "Tikh"], "Location", "eastoutside")
saveFigureAsEps(makeFileName(profileNumber, "model_parameter_error_comparison.eps"), fig)


%% Display Final Results

disp("Final Results")
disp(gyroTable)