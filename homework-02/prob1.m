%% Problem 1
%
% Exercise 1 from section 2.7

close all
clear
clc

warning off

% Add "Lib" to serach path
addpath(genpath(fullfile("..", "PEIP-master", "Lib")))

% Load examplar data
load(fullfile(pwd, "profile.mat"))

% Print Preamble
fprintf("Homework 2, Problem 1\n\n")


%% Part A

% Print preamble
fprintf("**************************\n")
fprintf("Part A\n")
fprintf("**************************\n\n")

% Number of obversations
m = length(x);

% From problem statement
sigma = 0.1;

% Build up the operator G
G = [ones(m, 1), x];

% Apply weighting
W = (1/sigma) .* eye(6);
Gw = W * G;
dw = W * t;

% Solve the least squares problem
m_L2 = (Gw.' * Gw) \ Gw.' * dw;
t_0 = m_L2(1);
s_2 = m_L2(2);

n = length(m_L2);
nu = m - n;

% Print results
fprintf("Weighted Least-Squares Solution\n")
fprintf("\tt_0: %8.6f \n", t_0)
fprintf("\ts_2: %8.6f \n", s_2)

% Compute residuals
r_L2 = dw - Gw*m_L2;

% Plot data, the fitted model, and the residuals
fig = figure("Name", "Part A - Data, Fitted Model, and Residuals");
tl = tiledlayout(2, 1, "Parent", fig);
title(tl, "Part A - Data, Fitted Model, and Residuals");

xs = 0 : 0.1 : 30;
fittedModel = t_0 + s_2 * xs;

ax = nexttile(1);
hold(ax, "on")
plot(xs, fittedModel, 'k')
plot(x, t, 'bo')
title("Data and Fitted Model")
xlabel("x [km]")
ylabel("t [sec]")
grid on
grid minor
legend(["Fitted Model", "Data"], "Location", "Best")

ax = nexttile(2);
hold(ax, "on")
plot(x, r_L2, 'rx')
title("Residuals")
xlabel("x [km]")
ylabel("residuals")
grid on
grid minor

linkaxes(tl.Children, 'x')

% Space out print outs
fprintf("\n\n")


%% Part B

% Print preamble
fprintf("**************************\n")
fprintf("Part B\n")
fprintf("**************************\n\n")

% Model covariance matrix
A = (Gw.' * Gw) \ Gw.';
C = A * A.';

% Print results
fprintf("Model Covariance Matrix\n")
disp(C)

% Model parameter correlation matrix
Rho = zeros(n, n);
for ii = 1 : n
    for jj = 1 : n
        Rho(ii,jj) = C(ii,jj) ./ sqrt(C(ii,ii) * C(jj,jj));
    end
end

% Print results
fprintf("Model Parameter Correlation Matrix\n")
disp(Rho)

% Comment on the results
fprintf("These model parameters have a high negative correlation!\n\n")

% Space out print outs
fprintf("\n")


%% Part C

% Print preamble
fprintf("**************************\n")
fprintf("Part C\n")
fprintf("**************************\n\n")

% Compute 95% confidence intervals
confidenceInterval95 = 1.96 .* sqrt(diag(C));
fprintf("95%% Confidence Intervals\n")
fprintf("\tt_0: %8.6f +/- %8.6f \n", m_L2(1), confidenceInterval95(1))
fprintf("\ts_2: %8.6f +/- %8.6f \n\n", m_L2(2), confidenceInterval95(2))

% Compute 95% confidence interval box
boxCorners = [...
    m_L2(1) - confidenceInterval95(1), m_L2(1) + confidenceInterval95(1); ...
    m_L2(2) - confidenceInterval95(2), m_L2(2) + confidenceInterval95(2)];
fprintf("95%% Confidence Interval Box\n")
disp(boxCorners)

% Compute ellipse parameters
theta=(0:.01:2*pi)';
r = zeros(length(theta), 2);
[u, lambda] = eig(inv(C));
deltaChiSq = chi2inv(0.95, nu);
delta = sqrt(deltaChiSq);
r(:,1) = ...
    (delta / sqrt(lambda(1,1))) * u(1,1) * cos(theta) + ...
    (delta / sqrt(lambda(2,2))) * u(1,2) * sin(theta);
r(:,2) = ...
    (delta / sqrt(lambda(1,1))) * u(2,1) * cos(theta) + ...
    (delta / sqrt(lambda(2,2))) * u(2,2) * sin(theta);

% Plot 95% Confidence Region and 95% Confidence Ellipse
fig = figure("Name", "Part C - 95% Confidence Region and Ellipse");
ax = gca;
hold(ax, "on")
f1 = fill(...
    [boxCorners(1,1), boxCorners(1,2), boxCorners(1,2), boxCorners(1,1)], ...
    [boxCorners(2,1), boxCorners(2,1), boxCorners(2,2), boxCorners(2,2)], ...
    "blue", ...
    "FaceAlpha", 0.4);
f2 = fill(...
    m_L2(1) + r(:,1), ...
    m_L2(2) + r(:,2), ...
    'm', ...
    'FaceAlpha', 0.4);
title("Part C - 95% Confidence Ellipse")
xlabel("t_0 [sec]")
ylabel("s_2 [km]")
grid on
grid minor
% axis equal
legend([f1, f2], ["95% Box", "95% Ellipse"], "Location", "eastoutside")

% Space out print outs
fprintf("\n")


%% Part D

% Print preamble
fprintf("**************************\n")
fprintf("Part D\n")
fprintf("**************************\n\n")

% Compute chi2 observed
chi2 = norm((dw - Gw*m_L2) ./ sigma).^2;
fprintf("Observed chi2: %8.6f \n", chi2)

% Compute p-value
p = 1 - chi2cdf(chi2, nu);
fprintf("p-value: %8.6f \n", p)

% Space out print outs
fprintf("\n\n")


%% Part E

% Print preamble
fprintf("**************************\n")
fprintf("Part E\n")
fprintf("**************************\n\n")

% Seed RNG
rng('default')

% Monte Carlo
numTrials = 1000;
mmc = zeros(numTrials, n);
chimc = zeros(numTrials, 1);
y_fit = G * m_L2;
for k = 1 : numTrials
    y_trial = y_fit + sigma.*randn(m,1);
    y_w_trial = y_trial ./ sigma;
    mmc(k,:) = (Gw \ y_w_trial).';
    chimc(k)= norm((G*mmc(k,:).' - y_trial) ./ sigma)^2;
end

% Model Monte Carlo Histograms
fig = figure("Name", "Part E - Model Parameter Monte Carlo");
tl = tiledlayout(2, 1, "Parent", fig);
title(tl, "Part E - Model Parameter MC Histograms")

ax = nexttile(1);
histogram(ax, mmc(:,1))
title("t_0 Histogram")
xlabel("Bins")
grid on
grid minor

ax = nexttile(2);
histogram(ax, mmc(:,2))
title("s_2 Histogram")
xlabel("Bins")
grid on
grid minor

% Chi^2 Monte Carlo Histograms
fig = figure("Name", "Part E - \chi^2 Monte Carlo Histogram");
histogram(gca, chimc)
title("Part E - \chi^2 Monte Carlo Histogram")
xlabel("Bins")
grid on
grid minor

% Space out print outs
fprintf("Something is wrong with these plots...\n")
fprintf("\n\n")


%% Part F

% Print preamble
fprintf("**************************\n")
fprintf("Part F\n")
fprintf("**************************\n\n")



% Space out print outs
fprintf("I can't evaluate this until Part E is fixed..\n")
fprintf("\n\n")


%% Part G

% Print preamble
fprintf("**************************\n")
fprintf("Part G\n")
fprintf("**************************\n\n")

% Solve for model parameters with IRLS with p = 1
m_L1 = irls(Gw, dw, eps, 1e-6, 1, 1e3);

% Print results
fprintf("IRLS Solution\n")
fprintf("\tt_0: %8.6f \n", m_L1(1))
fprintf("\ts_2: %8.6f \n", m_L1(2))

% Compute residuals
r_L1 = dw - Gw*m_L1;

% Plot data, the fitted model, and the residuals
fig = figure("Name", "Part G - Data, Fitted Model, and Residuals");
tl = tiledlayout(2, 1, "Parent", fig);
title(tl, "Part G - Data, Fitted Model, and Residuals");

xs = 0 : 0.1 : 30;
fittedModelL2 = m_L2(1) + m_L2(2)*xs;
fittedModelL1 = m_L1(1) + m_L1(2)*xs;

ax = nexttile(1);
hold(ax, "on")
plot(xs, fittedModelL2, 'k')
plot(xs, fittedModelL1, 'm')
plot(x, t, 'bo')
title("Data and Fitted Model")
xlabel("x [km]")
ylabel("t [sec]")
grid on
grid minor
legend(["L2 Fitted Model", "L1 Fitted Model", "Data"], "Location", "eastoutside")

ax = nexttile(2);
hold(ax, "on")
plot(x, r_L2, 'rx')
plot(x, r_L1, 'mx')
title("Residuals")
xlabel("x [km]")
ylabel("residuals")
grid on
grid minor
legend(["L2 Residuals", "L1 Residuals"], "Location", "eastoutside")

linkaxes(tl.Children, 'x')

% Space out print outs
fprintf("\n\n")


%% Part H

% Print preamble
fprintf("**************************\n")
fprintf("Part H\n")
fprintf("**************************\n\n")

% Monte Carlo
numTrials = 1000;
mmc_L1 = zeros(numTrials, n);
chimc_L1 = zeros(numTrials, 1);
y_fit = G * m_L1;
for k = 1 : numTrials
    y_trial = y_fit + sigma.*randn(m,1);
    y_w_trial = y_trial ./ sigma;
    mmc_L1(k,:) = irls(Gw, y_w_trial, eps, 1e-6, 1, 1e3);
    chimc_L1(k)= norm((G*mmc_L1(k,:).' - y_trial) ./ sigma)^2;
end

% Model Monte Carlo Histograms
fig = figure("Name", "Part H - L1 Model Parameter Monte Carlo");
tl = tiledlayout(2, 1, "Parent", fig);
title(tl, "Part H - L1 Model Parameter MC Histograms")

ax = nexttile(1);
histogram(ax, mmc_L1(:,1))
title("t_0 Histogram")
xlabel("Bins")
grid on
grid minor

ax = nexttile(2);
histogram(ax, mmc_L1(:,2))
title("s_2 Histogram")
xlabel("Bins")
grid on
grid minor

% Chi^2 Monte Carlo Histograms
fig = figure("Name", "Part H - L1 \chi^2 Monte Carlo Histogram");
histogram(gca, chimc_L1, 200)
title("Part H - L1 \chi^2 Monte Carlo Histogram")
xlabel("Bins")
grid on
grid minor

% Space out print outs
fprintf("\n\n")


%% Part I

% Print preamble
fprintf("**************************\n")
fprintf("Part I\n")
fprintf("**************************\n\n")



% Space out print outs
fprintf("\n\n")

