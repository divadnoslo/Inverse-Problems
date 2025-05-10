%% Task 11: Evaluate MP3 Accelerometer Experiment

% close all
% clear
% clc

% Save figures as *.eps
saveFigureAsEps = @(name, fig)(exportgraphics(fig, fullfile("..", "latex", "images", name)));

% Load Working File
load(fullfile(pwd, "working_file.mat"))

% Preamble
preamble = @(description)(sprintf("Motion Profile 3: %s", description));

% Make File Name
makeFileName = @(description)(sprintf("MP3_%s", description));


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
zx.YLim(1) = 0;
grid on
grid minor
saveFigureAsEps(makeFileName("accel_singular_values.eps"), fig)

% Estimate Calibration Parameters via Normal Equations
m_accel = inv(G.' * G) * G.' * d;
accelTable.L2Model = m_accel;

% Examine Error with True Model Parameters
m_accel_error = m_accel - m_accel_true;
accelTable.ModelError = m_accel_error;

% Model Covariance
C_accel = inv(G.' * G);

% Compute 95% Confidence Bound
conf95 = 1.96 * sqrt(diag(C_accel));
accelTable.Confidence95 = conf95;

% Model Correlation Matrix
Rho = computeModelCorrelationMatrix(C_accel);

% Plot Model Covariance Matrix
fig = figure("Name", "Accel Correlation matrix");
ax = gca;
colormap('gray')
imagesc(Rho)
clim([-1 1])
ax.XTick = 1 : n;
ax.YTick = 1 : n;
ax.XTickLabel = accelModelLabels;
ax.YTickLabel = accelModelLabels;
ax.YDir = "reverse";
axis equal
ax.XLim = [ax.XTick(1) - 0.5, ax.XTick(end) + 0.5];
ax.YLim = [ax.YTick(1) - 0.5, ax.YTick(end) + 0.5];
title(preamble("Accel Correlation Matrix"))
colorbar(ax, "eastoutside")
saveFigureAsEps(makeFileName("accel_correlation_matrix.eps"), fig)

% Display Accelermoeter Results
disp(preamble("Results"))
disp(accelTable)


%% Append to Working File

% save(workingFilePath, "w_b__i_b_true", "f_b__i_b_true", "w_b__i_b_meas", "f_b__i_b_meas", "-append")

