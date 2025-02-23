%% Problem 1
%
% Exercise 1 from section 2.7

close all 
clear
clc

warning off

% Add "Lib" to search path
addpath(genpath(fullfile("..", "PEIP-master", "Lib")))

% Save figures as *.eps
saveFigureAsEps = @(name, fig)(exportgraphics(fig, fullfile(pwd, "latex", "images", name)));

% Load examplar data
load(fullfile(pwd, "profile.mat"))

% Seed RNG
rng('default')

% Print preamble
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
W = (1/sigma) .* eye(m);
Gw = W * G;
dw = W * t;

% Solve the least squares problem
m_L2 = (Gw.' * Gw) \ Gw.' * dw;

n = length(m_L2);
nu = m - n;

% Print results
fprintf("Weighted Least-Squares Solution\n")
fprintf("\tt_0: %8.6f \n", m_L2(1))
fprintf("\ts_2: %8.6f \n", m_L2(2))

% Compute residuals
r_L2 = dw - Gw*m_L2;

% Plot data, the fitted model, and the residuals
fig = figure("Name", "Part A - Data, Fitted Model, and Residuals");
tl = tiledlayout(2, 1, "Parent", fig);
title(tl, "Part A - Data, Fitted Model, and Residuals");

xs = 0 : 0.1 : 30;
fittedModel = m_L2(1) + m_L2(2)*xs;

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
saveFigureAsEps("prob1_partA.eps", fig)

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
theta = (0 : 0.01 : 2*pi)';
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
saveFigureAsEps("prob1_partC.eps", fig)

% Space out print outs
fprintf("\n")


%% Part D

% Print preamble
fprintf("**************************\n")
fprintf("Part D\n")
fprintf("**************************\n\n")

% Compute chi2 observed
chi2 = sum((t - G*m_L2).^2 ./ sigma.^2);
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

% Monte Carlo
numTrials = 1000;
m_mc = zeros(numTrials, n);
chi_mc = zeros(numTrials, 1);
y_fit = G * m_L2;
for k = 1 : numTrials
    y_trial = y_fit + sigma.*randn(m,1);
    y_w_trial = y_trial ./ sigma;
    m_mc(k,:) = (Gw \ y_w_trial).';
    chi_mc(k)= norm((y_trial - G*m_mc(k,:).') ./ sigma^2);
end

% Model Monte Carlo Histograms
fig = figure("Name", "Part E - Model Parameter Monte Carlo");
tl = tiledlayout(2, 1, "Parent", fig);
title(tl, "Part E - Model Parameter MC Histograms")

ax = nexttile(1);
histogram(ax, m_mc(:,1))
title("t_0 Histogram")
xlabel("Bins")
grid on
grid minor

ax = nexttile(2);
histogram(ax, m_mc(:,2))
title("s_2 Histogram")
xlabel("Bins")
grid on
grid minor
saveFigureAsEps("prob1_partE_1.eps", fig)

% Chi^2 Monte Carlo Histograms
fig = figure("Name", "Part E - \chi^2 Monte Carlo Histogram");
histogram(gca, chi_mc)
title("Part E - \chi^2 Monte Carlo Histogram")
xlabel("Bins")
grid on
grid minor
saveFigureAsEps("prob1_partE_2.eps", fig)

% Space out print outs
fprintf("Plots only\n")
fprintf("\n\n")


%% Part F

% Print preamble
fprintf("**************************\n")
fprintf("Part F\n")
fprintf("**************************\n\n")

% Plot 95% Confidence Region and 95% Confidence Ellipse
fig = figure("Name", "Part F - Theoretical vs. MC");
ax = gca;
hold(ax, "on")
fill(...
    m_L2(1) + r(:,1), ...
    m_L2(2) + r(:,2), ...
    'm', ...
    'FaceAlpha', 0.4);
plot(m_mc(:,1), m_mc(:,2), 'k.')
title("Part F - Theoretical vs. Monte-Carlo Simulation")
xlabel("t_0 [sec]")
ylabel("s_2 [km]")
grid on
grid minor
% axis equal
legend(["95% Ellipse", "MC Realizations"], "Location", "best")
saveFigureAsEps("prob1_partF.eps", fig)

% Space out print outs
fprintf("Plots only\n")
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
plot(x, r_L1, 'rd')
title("Residuals")
xlabel("x [km]")
ylabel("residuals")
grid on
grid minor
legend(["L2 Residuals", "L1 Residuals"], "Location", "eastoutside")

linkaxes(tl.Children, 'x')
saveFigureAsEps("prob1_partG.eps", fig)

% Space out print outs
fprintf("\n\n")


%% Part H

% Print preamble
fprintf("**************************\n")
fprintf("Part H\n")
fprintf("**************************\n\n")

% Monte Carlo
numTrials = 1000;
m_mc_L1 = zeros(numTrials, n);
chi_mc_L1 = zeros(numTrials, 1);
y_fit = G * m_L1;
for k = 1 : numTrials
    y_trial = y_fit + sigma.*randn(m,1);
    y_w_trial = y_trial ./ sigma;
    m_mc_L1(k,:) = irls(Gw, y_w_trial, eps, 1e-6, 1, 1e3);
    chi_mc_L1(k)= norm((G*m_mc_L1(k,:).' - y_trial) ./ sigma)^2;
end

% Empirical estimate of the L1 regression covariance
A_L1 = m_mc_L1 - mean(m_mc_L1, 1);
C_L1 = (A_L1.' * A_L1) ./ numTrials;
fprintf("Empirical Estimate of Covariance for L1 Regression\n")
disp(C_L1)

% Resulting model parameter correlation matrix
Rho_L1 = zeros(n, n);
for ii = 1 : n
    for jj = 1 : n
        Rho_L1(ii,jj) = C_L1(ii,jj) ./ sqrt(C_L1(ii,ii) * C_L1(jj,jj));
    end
end
fprintf("Resulting Model Parameter Correlation Matrix\n")
disp(Rho_L1)

% Compute ellipse for comparison
theta = (0 : 0.01 : 2*pi)';
rr = zeros(length(theta), 2);
[u, lambda] = eig(inv(C_L1));
deltaChiSq = chi2inv(0.95, nu);
delta = sqrt(deltaChiSq);
rr(:,1) = ...
    (delta / sqrt(lambda(1,1))) * u(1,1) * cos(theta) + ...
    (delta / sqrt(lambda(2,2))) * u(1,2) * sin(theta);
rr(:,2) = ...
    (delta / sqrt(lambda(1,1))) * u(2,1) * cos(theta) + ...
    (delta / sqrt(lambda(2,2))) * u(2,2) * sin(theta);

% Compare 95% ellipse with MC realization
fig = figure("Name", "Part H - Theoretical vs. MC");
ax = gca;
hold(ax, "on")
fill(...
    m_L1(1) + rr(:,1), ...
    m_L1(2) + rr(:,2), ...
    'g', ...
    'FaceAlpha', 0.4);
plot(m_mc_L1(:,1), m_mc_L1(:,2), 'k.')
title("Part H - L1 Regression 95% Confience Interval Ellipse")
xlabel("t_0 [sec]")
ylabel("s_2 [km]")
grid on
grid minor
% axis equal
legend(["L1 95% Ellipse", "L1 MC Realizations"], "Location", "best")
saveFigureAsEps("prob1_partH.eps", fig)

% Space out print outs
fprintf("\n\n")


%% Part I

% Print preamble
fprintf("**************************\n")
fprintf("Part I\n")
fprintf("**************************\n\n")



% Space out print outs
fprintf("Response in LaTeX only\n")
fprintf("\n\n")

