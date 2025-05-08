close all
clear
clc

% Add to Path
addpath(genpath(fullfile("..", "PEIP-master", "Lib")))

% Save figures as *.eps
saveFigureAsEps = @(name, fig)(exportgraphics(fig, fullfile(pwd, "latex", "images", name)));

% Load Data
load(fullfile(pwd, "data", "be10.mat"))


%% Problem 1 - Exercise 10 in Section 11.6


%% Create Grid of Chi-Squared Values

% Givens
lambda = 4.998e-7;
P0 = 3.2;
mu = 0.0166;
z = depths;

% Function Handle
P = @(z)(P0 * exp(-mu*z));

% Points
d_ep = 1e-5;
d_T = 1000;
ep = 5e-6 : d_ep : 1e-3;
T = 500 : d_T : 199500; 

% Make Grid
[X, Y] = meshgrid(ep, T);

% Compute Chi-Squared over the Region
chi2 = zeros(size(X));
[rows, cols] = size(X);
for ii = 1 : rows
    for jj = 1 : cols
        chi2(ii,jj) = model2chi2([X(ii,jj); Y(ii,jj)]);
    end
end  

% Chi-Squared Surface
fig = figure("Name", "Chi-Squared Surface");
ax = gca;
hold(ax, "on")
sh = surf(X, Y, chi2);
title("\chi^2 Surface")
xlabel("\epsilon")
ylabel("T")
zlabel("\chi^2")
grid on
view(75, 35)
saveFigureAsEps("prob1_partA_chi_2_surface.eps", fig);


%% Part A

% Convert Chi^2 Values to Likelihoods
L = chi2likelihood(chi2);

% Likelihood Surface
fig = figure("Name", "Likelihood Surface");
ax = gca;
hold(ax, "on")
sh = surf(X, Y, L);
title("L(m|d) Surface")
xlabel("\epsilon")
ylabel("T")
zlabel("L(m|d)")
grid on
view(75, 35)
saveFigureAsEps("prob1_partA_likelihood_surface.eps", fig);

% Posterior Distribution
q = L ./ (sum(L, 'all') * d_ep * d_T);
ismembertol(sum(q, "all")  * d_ep * d_T, 1, 1e-3);

fig = figure("Name", "Posterior Distribution");
ax = gca;
hold(ax, "on")
sh = surf(X, Y, q);
title("Posterior Distribution q(m|d)")
xlabel("\epsilon")
ylabel("T")
zlabel("q(m|d)")
grid on
view(75, 35)
saveFigureAsEps("prob1_partA_posterior_distribution.eps", fig);


%% Part B

% Marginal Probabilities
p_ep = sum(q, 1) * d_T;
p_T = sum(q, 2) * d_ep;

% Ensure they are properly normalized
ismembertol(sum(p_ep * d_ep), 1, 1e-3);
ismembertol(sum(p_T * d_T), 1, 1e-3);

% Plot Marginal Distibutions
fig = figure("Name", "Marginal Distributions");
tl = tiledlayout(2, 1, "Parent", fig);
title(tl, "Marginal Distributions")
ax = nexttile(1);
area(ep, p_ep, 'FaceAlpha', 0.5)
title("p_\epsilon(\epsilon)")
xlabel("\epsilon")
ylabel("PDF")
grid on
grid minor
ax = nexttile(2);
area(T, p_T, 'FaceAlpha', 0.5)
title("p_T(T)")
xlabel("T")
ylabel("PDF")
grid on
grid minor
saveFigureAsEps("prob1_partB_marginal_distributions.eps", fig);


%% Part C

% Prior Distribution
prior = NormalDistribution(0.0005, 0.0002^2);

% Plot Prior
fig = figure("Name", "Prior Distribution");
area(ep, prior.probabilityDensityFunction(ep), 'FaceAlpha', 0.5)
title("Prior Distribution")
xlabel("\epsilon")
ylabel("p(\epsilon)")
grid on
grid minor
saveFigureAsEps("prob1_partC_prior_distribution.eps", fig);

% Apply Prior to the Likelihood
qP = zeros(size(X));
[rows, cols] = size(X);
for ii = 1 : rows
    for jj = 1 : cols
        qP(ii,jj) = prior.probabilityDensityFunction(X(ii,jj)) * L(ii,jj);
    end
end  

% Normalize
qN = qP ./ (sum(qP, 'all') * d_ep * d_T);
ismembertol(sum(qN, "all")  * d_ep * d_T, 1, 1e-3);

fig = figure("Name", "Posterior Distribution with Prior Applied");
ax = gca;
hold(ax, "on")
sh = surf(X, Y, qN);
title("Posterior Distribution q(m|d) with Prior Applied")
xlabel("\epsilon")
ylabel("T")
zlabel("q(m|d)")
grid on
view(75, 35)
saveFigureAsEps("prob1_partC_posterior_distribution.eps", fig);

% Marginal Probabilities
p_ep_N = sum(qN, 1) * d_T;
p_T_N = sum(qN, 2) * d_ep;

% Ensure they are properly normalized
ismembertol(sum(p_ep_N * d_ep), 1, 1e-3);
ismembertol(sum(p_T_N * d_T), 1, 1e-3);

% Plot Marginal Distibutions
fig = figure("Name", "Marginal Distributions");
tl = tiledlayout(2, 1, "Parent", fig);
title(tl, "Marginal Distributions with Prior Applied")
ax = nexttile(1);
area(ep, p_ep_N, 'FaceAlpha', 0.5)
title("p_\epsilon(\epsilon)")
xlabel("\epsilon")
ylabel("PDF")
grid on
grid minor
ax = nexttile(2);
area(T, p_T_N, 'FaceAlpha', 0.5)
title("p_T(T)")
xlabel("T")
ylabel("PDF")
grid on
grid minor
saveFigureAsEps("prob1_partC_marginal_distributions.eps", fig);