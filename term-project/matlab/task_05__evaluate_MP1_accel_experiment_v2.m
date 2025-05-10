%% Task 04: Evaluate Single-Axis Motion Accelerometer Experiment V2

% close all
% clear
% clc

% Add "Lib" to Path
addpath(genpath(fullfile("..", "..", "PEIP-master", "Lib")))

% Save figures as *.eps
saveFigureAsEps = @(name, fig)(exportgraphics(fig, fullfile("..", "latex", "images", name)));

% Load Working File
load(fullfile(pwd, "working_file.mat"))

% Preamble
preamble = @(description)(sprintf("Motion Profile 1: %s", description));

% Make File Name
makeFileName = @(description)(sprintf("MP1_%s", description));


%% Evaluate Accelerometer Calibration Parameters

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
title(preamble("Accelerometer Singular Values"))
xlabel("Singular Value Index")
ylabel("s_i")
grid on
grid minor
saveFigureAsEps(makeFileName("accel_singular_values_v2.eps"), fig)

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
accelTable.SvdModel = m_svd;

% Model Error
accelTable.SvdModelError = m_svd - m_accel_true;

% Model Null Space
disp("Model Null Space")
disp(V0)


%% Tikohonov Regularization

s = diag(S);

% L-Curve
[rho, eta, reg_param] = l_curve_tikh_svd(U, s, d, 1000, 1e-12, 1e12);
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
title(preamble("Zeroth-Order Tikhonov L-Curve"))
xlabel("||Gm - d||_2")
ylabel("||m||_2")
grid on
grid minor
saveFigureAsEps(makeFileName("accel_0th_order_tik_l_curve.eps"), fig)

% Compute Zeroth-Order Tikhonov Solution
m_tikh = (G.'*G + (alpha_tikh^2)*eye(n)) \ G.' * d;
accelTable.TikhModel = m_tikh;

% Model Error
accelTable.TikhModelError = m_tikh - m_accel_true;


%% Display Final Results

disp("Final Results")
disp(accelTable)

