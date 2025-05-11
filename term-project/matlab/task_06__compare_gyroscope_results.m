%% Task 06: Compare Gyroscope Results

% close all
% clear
% clc

% Load Working File
load(fullfile(pwd, "working_file.mat"))

% Helper Functions
saveFigureAsEps = @(name, fig)(exportgraphics(fig, fullfile("..", "latex", "images", name)));


%% Unpack Results

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
mp1Color = "#FF0000";
mp2Color = "#00FF00";
mp3Color = "#0000FF";


%% Model Parameter Results

% Unpack Model Parameters
gyroBiasResults = [...
    m_gyro_true(bIndexes).'; ...
    gyroResults(1).results.L2Model(bIndexes).'; ...
    gyroResults(2).results.L2Model(bIndexes).'; ...
    gyroResults(3).results.L2Model(bIndexes).'].';
gyroScaleFactorResults = [...
    m_gyro_true(sfIndexes).'; ...
    gyroResults(1).results.L2Model(sfIndexes).'; ...
    gyroResults(2).results.L2Model(sfIndexes).'; ...
    gyroResults(3).results.L2Model(sfIndexes).'].';
gyroMisalignmentResults = [...
    m_gyro_true(mIndexes).'; ...
    gyroResults(1).results.L2Model(mIndexes).'; ...
    gyroResults(2).results.L2Model(mIndexes).'; ...
    gyroResults(3).results.L2Model(mIndexes).'];

% Model Results Bar Chart
fig = figure("Name", "Gyro Model Bar Chart");
tl = tiledlayout(3, 1, "Parent", fig);
title(tl, "Gyroscope Calibration Parameter Estimation Results")
ax = nexttile(1);
b = bar(gyroBiases, 180/pi * gyroBiasResults);
title("Bias Estimation Comparison")
ylabel("[deg/sec]")
b(1).FaceColor = trueColor;
b(2).FaceColor = mp1Color;
b(3).FaceColor = mp2Color;
b(4).FaceColor = mp3Color;
grid on
grid minor
legend(["Truth", "MP1", "MP2", "MP3"], "Location", "eastoutside")
ax = nexttile(2);
b = bar(gyroScaleFactors, 1e6 * gyroScaleFactorResults);
title("Scale Factor Estimation Comparison")
ylabel("[ppm]")
b(1).FaceColor = trueColor;
b(2).FaceColor = mp1Color;
b(3).FaceColor = mp2Color;
b(4).FaceColor = mp3Color;
grid on
grid minor
legend(["Truth", "MP1", "MP2", "MP3"], "Location", "eastoutside")
ax = nexttile(3);
b = bar(gyroMisalignments, 1e3 * gyroMisalignmentResults);
title("Misalignment Estimation Comparison")
ylabel("[m-rad]")
b(1).FaceColor = trueColor;
b(2).FaceColor = mp1Color;
b(3).FaceColor = mp2Color;
b(4).FaceColor = mp3Color;
grid on
grid minor
legend(["Truth", "MP1", "MP2", "MP3"], "Location", "eastoutside")
saveFigureAsEps("gyro_parameter_comparison.eps", fig)


%% Model Errors

% Unpack Model Errors
gyroBiasErrors = [...
    gyroResults(1).results.ModelError(bIndexes).'; ...
    gyroResults(2).results.ModelError(bIndexes).'; ...
    gyroResults(3).results.ModelError(bIndexes).'].';
gyroScaleFactorErrors = [...
    gyroResults(1).results.ModelError(sfIndexes).'; ...
    gyroResults(2).results.ModelError(sfIndexes).'; ...
    gyroResults(3).results.ModelError(sfIndexes).'].';
gyroMisalignmentErrors = [...
    gyroResults(1).results.ModelError(mIndexes).'; ...
    gyroResults(2).results.ModelError(mIndexes).'; ...
    gyroResults(3).results.ModelError(mIndexes).'];

% Model Error Bar Chart
fig = figure("Name", "Gyro Error Bar Chart");
tl = tiledlayout(3, 1, "Parent", fig);
title(tl, "Gyroscope Calibration Parameter Absolute Error")
ax = nexttile(1);
b = bar(gyroBiases, 180/pi * 1e3 * abs(gyroBiasErrors));
title("Bias Estimation Absolute Error")
ylabel("[milli-deg/sec]")
b(1).FaceColor = mp1Color;
b(2).FaceColor = mp2Color;
b(3).FaceColor = mp3Color;
grid on
grid minor
legend(["MP1", "MP2", "MP3"], "Location", "eastoutside")
ax = nexttile(2);
b = bar(gyroScaleFactors, 1e6 * abs(gyroScaleFactorErrors));
title("Scale Factor Estimation Absolute Error")
ylabel("[ppm]")
b(1).FaceColor = mp1Color;
b(2).FaceColor = mp2Color;
b(3).FaceColor = mp3Color;
grid on
grid minor
legend(["MP1", "MP2", "MP3"], "Location", "eastoutside")
ax = nexttile(3);
b = bar(gyroMisalignments, 1e6 * abs(gyroMisalignmentErrors));
title("Misalignment Estimation Absolute Error")
ylabel("[micro-rad]")
b(1).FaceColor = mp1Color;
b(2).FaceColor = mp2Color;
b(3).FaceColor = mp3Color;
grid on
grid minor
legend(["MP1", "MP2", "MP3"], "Location", "eastoutside")
saveFigureAsEps("gyro_parameter_error_comparison.eps", fig)


%% Covariances

% Unpack Model Covariances
gyroBiasCovariances = [...
    gyroResults(1).results.Covariance(bIndexes).'; ...
    gyroResults(2).results.Covariance(bIndexes).'; ...
    gyroResults(3).results.Covariance(bIndexes).'].';
gyroScaleFactorCovariances = [...
    gyroResults(1).results.Covariance(sfIndexes).'; ...
    gyroResults(2).results.Covariance(sfIndexes).'; ...
    gyroResults(3).results.Covariance(sfIndexes).'].';
gyroMisalignmentCovariances = [...
    gyroResults(1).results.Covariance(mIndexes).'; ...
    gyroResults(2).results.Covariance(mIndexes).'; ...
    gyroResults(3).results.Covariance(mIndexes).'];

% Model Covariances Bar Chart
fig = figure("Name", "Gyro Error Bar Chart");
tl = tiledlayout(3, 1, "Parent", fig);
title(tl, "Gyroscope Calibration Parameter Covariances")
ax = nexttile(1);
b = bar(gyroBiases, 180/pi * gyroBiasCovariances);
title("Bias Estimation Covariances")
ylabel("[milli-deg/sec]")
b(1).FaceColor = mp1Color;
b(2).FaceColor = mp2Color;
b(3).FaceColor = mp3Color;
grid on
grid minor
legend(["MP1", "MP2", "MP3"], "Location", "eastoutside")
ax = nexttile(2);
b = bar(gyroScaleFactors, 1e6 * gyroScaleFactorCovariances);
title("Scale Factor Estimation Covariances")
ylabel("[ppm]")
ax.YScale = "log";
% ax.YLim(1) = 1e1;
b(1).FaceColor = mp1Color;
b(2).FaceColor = mp2Color;
b(3).FaceColor = mp3Color;
grid on
grid minor
legend(["MP1", "MP2", "MP3"], "Location", "eastoutside")
ax = nexttile(3);
b = bar(gyroMisalignments, 1e3 * gyroMisalignmentCovariances);
title("Misalignment Estimation Covariances")
ylabel("[milli-rad]")
ax.YScale = "log";
% ax.YLim(1) = 1e-2;
b(1).FaceColor = mp1Color;
b(2).FaceColor = mp2Color;
b(3).FaceColor = mp3Color;
grid on
grid minor
legend(["MP1", "MP2", "MP3"], "Location", "eastoutside")
saveFigureAsEps("gyro_parameter_covariance_comparison.eps", fig)


%% Confidence Intervals

% Unpack Model Confidence Intervals
gyroBiasConf95 = [...
    gyroResults(1).results.Confidence95(bIndexes).'; ...
    gyroResults(2).results.Confidence95(bIndexes).'; ...
    gyroResults(3).results.Confidence95(bIndexes).'].';
gyroScaleFactorConf95 = [...
    gyroResults(1).results.Confidence95(sfIndexes).'; ...
    gyroResults(2).results.Confidence95(sfIndexes).'; ...
    gyroResults(3).results.Confidence95(sfIndexes).'].';
gyroMisalignmentConf95 = [...
    gyroResults(1).results.Confidence95(mIndexes).'; ...
    gyroResults(2).results.Confidence95(mIndexes).'; ...
    gyroResults(3).results.Confidence95(mIndexes).'];

% Model Confidence Intervals Bar Chart
fig = figure("Name", "Gyro Error Bar Chart");
tl = tiledlayout(3, 1, "Parent", fig);
title(tl, "Gyroscope Calibration Parameter 95% Confidence Interval ")
ax = nexttile(1);
b = bar(gyroBiases, 180/pi * gyroBiasConf95);
title("Bias")
ylabel("[milli-deg/sec]")
b(1).FaceColor = mp1Color;
b(2).FaceColor = mp2Color;
b(3).FaceColor = mp3Color;
grid on
grid minor
legend(["MP1", "MP2", "MP3"], "Location", "eastoutside")
ax = nexttile(2);
b = bar(gyroScaleFactors, 1e6 * gyroScaleFactorConf95);
title("Scale Factor")
ylabel("[ppm]")
b(1).FaceColor = mp1Color;
b(2).FaceColor = mp2Color;
b(3).FaceColor = mp3Color;
grid on
grid minor
legend(["MP1", "MP2", "MP3"], "Location", "eastoutside")
ax = nexttile(3);
b = bar(gyroMisalignments, 1e3 * gyroMisalignmentConf95);
title("Misalignment")
ylabel("[milli-rad]")
b(1).FaceColor = mp1Color;
b(2).FaceColor = mp2Color;
b(3).FaceColor = mp3Color;
grid on
grid minor
legend(["MP1", "MP2", "MP3"], "Location", "eastoutside")
saveFigureAsEps("gyro_parameter_conf95_comparison.eps", fig)