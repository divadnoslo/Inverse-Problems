%% Task 03: Evaluate Motion Profile 1 -  Gyroscope Experiment

% close all
% clear
% clc

% Save figures as *.eps
saveFigureAsEps = @(name, fig)(exportgraphics(fig, fullfile("..", "latex", "images", name)));

% Load Working File
load(fullfile(pwd, "working_file.mat"))

% Preamble
preamble = @(description)(sprintf("Motion Profile 1: %s", description));

% Make File Name
makeFileName = @(description)(sprintf("MP1_%s", description));


%% Evaluate Gyropscope Calibration Parameters

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
% W = gyroSigma * eye(m);
% G = W * G;
% d = W * d;

% Examine Singular Values
[U, S, V] = svd(G, 'econ'); 
fig = figure("Name", "SVD Singular Values");
ax = gca;
hold(ax, "on")
plot(1:n, diag(S), 'bo')
title(preamble("Gyroscope Singular Values"))
xlabel("Singular Value Index")
ylabel("s_i")
ax.YLim(1) = 0;
ax.YScale = "log";
grid on
grid minor
saveFigureAsEps(makeFileName("gyro_singular_values.eps"), fig)

% Estimate Calibration Parameters via Normal Equations
m_gyro = inv(G.' * G) * G.' * d;
gyroTable.L2Model = m_gyro;

% Examine Error with True Model Parameters
m_gyro_error = m_gyro - m_gyro_true;
gyroTable.ModelError = m_gyro_error;

% Model Covariance
C_gyro = inv(G.' * G);

% Compute 95% Confidence Bound
conf95 = 1.96 * sqrt(diag(C_gyro));
gyroTable.Confidence95 = conf95;

% Model Correlation Matrix
Rho = computeModelCorrelationMatrix(C_gyro);

% Plot Model Covariance Matrix
fig = figure("Name", "Gyro Correlation matrix");
ax = gca;
colormap('gray')
imagesc(Rho)
clim([-1 1])
ax.XTick = 1 : n;
ax.YTick = 1 : n;
ax.XTickLabel = gyroModelLabels;
ax.YTickLabel = gyroModelLabels;
ax.YDir = "reverse";
axis equal
ax.XLim = [ax.XTick(1) - 0.5, ax.XTick(end) + 0.5];
ax.YLim = [ax.YTick(1) - 0.5, ax.YTick(end) + 0.5];
title(preamble("Gyro Correlation Matrix"))
colorbar(ax, "eastoutside")
saveFigureAsEps(makeFileName("gyro_correlation_matrix.eps"), fig)

% Display Gyro Results
disp(preamble("Results"))
disp(gyroTable)


%% Append to Working File

% save(workingFilePath, "gyroTable", "-append")

