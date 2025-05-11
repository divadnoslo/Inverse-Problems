%% Task 02: Run Traditional Motion Profile

% close all
% clear
% clc

% Load Working File
load(fullfile(pwd, "working_file.mat"))

% Profile Number
profileNumber = 0;

% Helper Functions
saveFigureAsEps = @(name, fig)(exportgraphics(fig, fullfile("..", "latex", "images", name)));
preamble = @(profileNumber, description)(sprintf("Motion Profile %s: %s", num2str(profileNumber), description));
makeFileName = @(profileNumber, description)(sprintf("MP%s_%s", num2str(profileNumber), description));


%% Simulate Traditional Motion Profile

imuCalManager = ImuCalibrationManager();

calData = imuCalManager.createCalibrationDataset(imu);
% imuCalManager.plotCalibrationDataset(calData);
processedImuModel = imuCalManager.processCalibrationDataset(calData);


%% Construct Model Operators

[G_g, d_g] = dataset2gyroInverseProblem(calData);
[G_a, d_a] = dataset2gyroInverseProblem(calData);


%% Evaluate Gyroscope Calibration

G = G_g;
d = d_g;
[m, n] = size(G);

% Examine Singular Values
[U, S, V] = svd(G, 'econ'); 
fig = figure("Name", "SVD Singular Values");
ax = gca;
hold(ax, "on")
plot(1:n, diag(S), 'bo')
title(preamble(profileNumber, "Gyroscope Singular Values"))
xlabel("Singular Value Index")
ylabel("s_i")
ax.XLim(1) = 1;
ax.XLim(2) = n;
ax.YLim(1) = 1e-1;
ax.YScale = "log";
grid on
grid minor
saveFigureAsEps(makeFileName(profileNumber, "gyro_singular_values.eps"), fig)

% Check if G' * G is invertable
if ismembertol(det(G.' * G), 0, 1e-6)
    fprintf("WARNING: (G_g *.' G_g) is not invertable!\n")
end


%% Evaluate Accelerometer Calibration

G = G_a;
d = d_a;
[m, n] = size(G);

% Examine Singular Values
[U, S, V] = svd(G, 'econ'); 
fig = figure("Name", "SVD Singular Values");
ax = gca;
hold(ax, "on")
plot(1:n, diag(S), 'bo')
title(preamble(profileNumber, "Accelerometer Singular Values"))
xlabel("Singular Value Index")
ylabel("s_i")
ax.XLim(1) = 1;
ax.XLim(2) = n;
ax.YLim(1) = 1e-1;
ax.YScale = "log";
grid on
grid minor
saveFigureAsEps(makeFileName(profileNumber, "accel_singular_values.eps"), fig)

% Check if G' * G is invertable
if ismembertol(det(G.' * G), 0, 1e-6)
    fprintf("WARNING: (G_a *.' G_a) is not invertable!\n\n")
end


%% Make Space

fprintf("\n\n")