%% Task 06: Compare Accelerometer Results

% close all
% clear
% clc

% Load Working File
load(fullfile(pwd, "working_file.mat"))

% Helper Functions
saveFigureAsEps = @(name, fig)(exportgraphics(fig, fullfile("..", "latex", "images", name)));


%% Unpack Results

% Put Into One Struct
accelResults(1) = mp1Accel;
accelResults(2) = mp2Accel;
accelResults(3) = mp3Accel;

% Indexes
bIndexes = [1, 5, 9];
sfIndexes = [2, 7, 12];
mIndexes = [3, 4, 6, 8, 10, 11];

% Create Catagories
accelBiases = categorical(accelModelLabels(bIndexes));
accelScaleFactors = categorical(accelModelLabels(sfIndexes));
accelMisalignments = categorical(accelModelLabels(mIndexes));

% Set Colors
trueColor = "#000000";
mp1Color = "#FF0000";
mp2Color = "#00FF00";
mp3Color = "#0000FF";


%% Model Parameter Results

% Unpack Model Parameters
accelBiasResults = [...
    m_accel_true(bIndexes).'; ...
    accelResults(1).results.L2Model(bIndexes).'; ...
    accelResults(2).results.L2Model(bIndexes).'; ...
    accelResults(3).results.L2Model(bIndexes).'].';
accelScaleFactorResults = [...
    m_accel_true(sfIndexes).'; ...
    accelResults(1).results.L2Model(sfIndexes).'; ...
    accelResults(2).results.L2Model(sfIndexes).'; ...
    accelResults(3).results.L2Model(sfIndexes).'].';
accelMisalignmentResults = [...
    m_accel_true(mIndexes).'; ...
    accelResults(1).results.L2Model(mIndexes).'; ...
    accelResults(2).results.L2Model(mIndexes).'; ...
    accelResults(3).results.L2Model(mIndexes).'];

% Model Results Bar Chart
fig = figure("Name", "Accel Model Bar Chart");
tl = tiledlayout(3, 1, "Parent", fig);
title(tl, "Accelerometer Calibration Parameter Estimation Results")
ax = nexttile(1);
b = bar(accelBiases, accelBiasResults);
title("Bias Estimation Comparison")
ylabel("[m/s^2]")
b(1).FaceColor = trueColor;
b(2).FaceColor = mp1Color;
b(3).FaceColor = mp2Color;
b(4).FaceColor = mp3Color;
grid on
grid minor
legend(["Truth", "MP1", "MP2", "MP3"], "Location", "eastoutside")
ax = nexttile(2);
b = bar(accelScaleFactors, 1e6 * accelScaleFactorResults);
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
b = bar(accelMisalignments, 1e3 * accelMisalignmentResults);
title("Misalignment Estimation Comparison")
ylabel("[m-rad]")
b(1).FaceColor = trueColor;
b(2).FaceColor = mp1Color;
b(3).FaceColor = mp2Color;
b(4).FaceColor = mp3Color;
grid on
grid minor
legend(["Truth", "MP1", "MP2", "MP3"], "Location", "eastoutside")
saveFigureAsEps("accel_parameter_comparison.eps", fig)


%% Chi^2 Comparison

% Unpack Chi^2 Values
chi2 = zeros(3, 1);
for k = 1 : 3
    chi2(k) = norm(accelResults(k).results.ModelError)^2;
end

% Create Motion Profile Categories
mps = categorical({'MP1', 'MP2', 'MP3'});

% Model Error Bar Chart
fig = figure("Name", "Accel Residuals Chart");
ax = gca;
b = bar(mps, chi2, 'FaceColor', 'flat');
b.CData(1,:) = [1 0 0];
b.CData(2,:) = [0 1 0];
b.CData(3,:) = [0 0 1];
title("Accelerometer \chi^2 Value Comparison")
ylabel("\chi^2")
ax.YScale = "log";
grid on
grid minor
saveFigureAsEps("accel_chi2_comparison.eps", fig)


%% Model Errors

% Unpack Model Errors
accelBiasErrors = [...
    accelResults(1).results.ModelError(bIndexes).'; ...
    accelResults(2).results.ModelError(bIndexes).'; ...
    accelResults(3).results.ModelError(bIndexes).'].';
accelScaleFactorErrors = [...
    accelResults(1).results.ModelError(sfIndexes).'; ...
    accelResults(2).results.ModelError(sfIndexes).'; ...
    accelResults(3).results.ModelError(sfIndexes).'].';
accelMisalignmentErrors = [...
    accelResults(1).results.ModelError(mIndexes).'; ...
    accelResults(2).results.ModelError(mIndexes).'; ...
    accelResults(3).results.ModelError(mIndexes).'];

% Model Error Bar Chart
fig = figure("Name", "Accel Error Bar Chart");
tl = tiledlayout(3, 1, "Parent", fig);
title(tl, "Accelerometer Calibration Parameter Absolute Error")
ax = nexttile(1);
b = bar(accelBiases, abs(accelBiasErrors));
title("Bias Estimation Absolute Error")
ylabel("[m/s^2]")
ax.YScale = "log";
ax.YLim(1) = 1e-7;
b(1).FaceColor = mp1Color;
b(2).FaceColor = mp2Color;
b(3).FaceColor = mp3Color;
grid on
grid minor
legend(["MP1", "MP2", "MP3"], "Location", "eastoutside")
ax = nexttile(2);
b = bar(accelScaleFactors, 1e6 * abs(accelScaleFactorErrors));
title("Scale Factor Estimation Absolute Error")
ylabel("[ppm]")
ax.YScale = "log";
b(1).FaceColor = mp1Color;
b(2).FaceColor = mp2Color;
b(3).FaceColor = mp3Color;
grid on
grid minor
legend(["MP1", "MP2", "MP3"], "Location", "eastoutside")
ax = nexttile(3);
ax.YScale = "log";
b = bar(accelMisalignments, 1e6 * abs(accelMisalignmentErrors));
title("Misalignment Estimation Absolute Error")
ylabel("[micro-rad]")
ax.YScale = "log";
ax.YLim(1) = 1e-1;
b(1).FaceColor = mp1Color;
b(2).FaceColor = mp2Color;
b(3).FaceColor = mp3Color;
grid on
grid minor
legend(["MP1", "MP2", "MP3"], "Location", "eastoutside")
saveFigureAsEps("accel_parameter_error_comparison.eps", fig)


%% Covariances

% Unpack Model Covariances
accelBiasCovariances = [...
    accelResults(1).results.Covariance(bIndexes).'; ...
    accelResults(2).results.Covariance(bIndexes).'; ...
    accelResults(3).results.Covariance(bIndexes).'].';
accelScaleFactorCovariances = [...
    accelResults(1).results.Covariance(sfIndexes).'; ...
    accelResults(2).results.Covariance(sfIndexes).'; ...
    accelResults(3).results.Covariance(sfIndexes).'].';
accelMisalignmentCovariances = [...
    accelResults(1).results.Covariance(mIndexes).'; ...
    accelResults(2).results.Covariance(mIndexes).'; ...
    accelResults(3).results.Covariance(mIndexes).'];

% Model Covariances Bar Chart
fig = figure("Name", "Accel Error Bar Chart");
tl = tiledlayout(3, 1, "Parent", fig);
title(tl, "Accelerometer Calibration Parameter Covariances")
ax = nexttile(1);
b = bar(accelBiases, accelBiasCovariances);
title("Bias Estimation Covariances")
ylabel("[m/s^2]")
ax.YScale = "log";
ax.YLim(1) = 1e-4;
b(1).FaceColor = mp1Color;
b(2).FaceColor = mp2Color;
b(3).FaceColor = mp3Color;
grid on
grid minor
legend(["MP1", "MP2", "MP3"], "Location", "eastoutside")
ax = nexttile(2);
ax.YScale = "log";
b = bar(accelScaleFactors, 1e6 * accelScaleFactorCovariances);
title("Scale Factor Estimation Covariances")
ylabel("[ppm]")
ax.YScale = "log";
ax.YLim(1) = 1e1;
b(1).FaceColor = mp1Color;
b(2).FaceColor = mp2Color;
b(3).FaceColor = mp3Color;
grid on
grid minor
legend(["MP1", "MP2", "MP3"], "Location", "eastoutside")
ax = nexttile(3);
ax.YScale = "log";
b = bar(accelMisalignments, 1e3 * accelMisalignmentCovariances);
title("Misalignment Estimation Covariances")
ylabel("[milli-rad]")
ax.YScale = "log";
ax.YLim(1) = 1e-3;
b(1).FaceColor = mp1Color;
b(2).FaceColor = mp2Color;
b(3).FaceColor = mp3Color;
grid on
grid minor
legend(["MP1", "MP2", "MP3"], "Location", "eastoutside")
saveFigureAsEps("accel_parameter_covariance_comparison.eps", fig)


%% Confidence Intervals

% Unpack Model Confidence Intervals
accelBiasConf95 = [...
    accelResults(1).results.Confidence95(bIndexes).'; ...
    accelResults(2).results.Confidence95(bIndexes).'; ...
    accelResults(3).results.Confidence95(bIndexes).'].';
accelScaleFactorConf95 = [...
    accelResults(1).results.Confidence95(sfIndexes).'; ...
    accelResults(2).results.Confidence95(sfIndexes).'; ...
    accelResults(3).results.Confidence95(sfIndexes).'].';
accelMisalignmentConf95 = [...
    accelResults(1).results.Confidence95(mIndexes).'; ...
    accelResults(2).results.Confidence95(mIndexes).'; ...
    accelResults(3).results.Confidence95(mIndexes).'];

% Model Confidence Intervals Bar Chart
fig = figure("Name", "Accel Error Bar Chart");
tl = tiledlayout(3, 1, "Parent", fig);
title(tl, "Accelerometer Calibration Parameter 95% Confidence Interval ")
ax = nexttile(1);
b = bar(accelBiases, accelBiasConf95);
title("Bias")
ylabel("[m/s^2]")
b(1).FaceColor = mp1Color;
b(2).FaceColor = mp2Color;
b(3).FaceColor = mp3Color;
grid on
grid minor
legend(["MP1", "MP2", "MP3"], "Location", "eastoutside")
ax = nexttile(2);
b = bar(accelScaleFactors, 1e6 * accelScaleFactorConf95);
title("Scale Factor")
ylabel("[ppm]")
b(1).FaceColor = mp1Color;
b(2).FaceColor = mp2Color;
b(3).FaceColor = mp3Color;
grid on
grid minor
legend(["MP1", "MP2", "MP3"], "Location", "eastoutside")
ax = nexttile(3);
b = bar(accelMisalignments, 1e3 * accelMisalignmentConf95);
title("Misalignment")
ylabel("[milli-rad]")
b(1).FaceColor = mp1Color;
b(2).FaceColor = mp2Color;
b(3).FaceColor = mp3Color;
grid on
grid minor
legend(["MP1", "MP2", "MP3"], "Location", "eastoutside")
saveFigureAsEps("accel_parameter_conf95_comparison.eps", fig)