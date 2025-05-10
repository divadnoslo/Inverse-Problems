close all
clear
clc

tic

% Add to Path
addpath(genpath(fullfile("..", "PEIP-master", "Lib")))

% Save figures as *.eps
saveFigureAsEps = @(name, fig)(exportgraphics(fig, fullfile(pwd, "latex", "images", name)));

% Load Data
load(fullfile(pwd, "data", "be10.mat"))


%% Problem 2 - Exercise 11 in Section 11.6

% Points
d_ep = 1e-5;
d_T = 1000;
ep = 5e-6 : d_ep : 1e-3;
T = 500 : d_T : 199500; 

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

%% MCMC Settings

global stepsize

N = 1e4;
stepsize = 0.5 * [d_ep; d_T];


%% Part A

m0 = [...
    ep(50); ...
    T(100)];
[mout,mMAP,pacc]=mcmc('logprior1','loglikelihood','generate','logproposal',m0,N);

disp("Part A")
disp(mMAP)

% MAP Solution 1
fig = figure("Name", "M-MAP Solution 1");
ax = gca;
hold(ax, "on")
plot(mout(1,:), mout(2,:), 'k.')
plot(m0(1), m0(2), 'b.', 'MarkerSize', 20)
plot(mMAP(1), mMAP(2), 'r.', 'MarkerSize', 20)
title("MCMC \rightarrow m_M_A_P: Uniform Prior")
xlabel("\epsilon")
ylabel("T")
xlim([ep(1), max(mout(1,:) + d_ep)])
ylim([T(1) T(end)])
grid on
legend(["Candidate Model", "Initial Model", "m_M_A_P"], 'Location', 'eastoutside')
saveFigureAsEps("prob2_uniform_prior.eps", fig);


%% Part B

% Plot Marginal Probabilities
fig = figure("Name", "Marginal Probabilities with Uniform Prior");
tl = tiledlayout(2, 1, "Parent", fig);
title(tl, "Marginal Probabilities with Uniform Prior")
ax = nexttile(1);
histogram(mout(1,:), "Normalization", "pdf")
title("p_\epsilon(\epsilon)")
xlabel("\epsilon")
ylabel("PDF")
xlim([ep(1) ep(end)])
grid on
grid minor
ax = nexttile(2);
histogram(mout(2,:), "Normalization", "pdf")
title("p_T(T)")
xlabel("T")
ylabel("PDF")
xlim([T(1) T(end)])
grid on
grid minor
saveFigureAsEps("prob2_uniform_prior_marginal_prob.eps", fig);


%% Part C

[mout,mMAP,pacc]=mcmc('logprior2','loglikelihood','generate','logproposal',m0,N);

disp("Part B")
disp(mMAP)

% MAP Solution 2
fig = figure("Name", "M-MAP Solution 2");
ax = gca;
hold(ax, "on")
plot(mout(1,:), mout(2,:), 'k.')
plot(m0(1), m0(2), 'b.', 'MarkerSize', 20)
plot(mMAP(1), mMAP(2), 'r.', 'MarkerSize', 20)
title("MCMC \rightarrow m_M_A_P: Normal Prior")
xlabel("\epsilon")
ylabel("T")
xlim([ep(1) ep(end)])
ylim([T(1) T(end)])
grid on
legend(["Candidate Model", "Initial Model", "m_M_A_P"], 'Location', 'eastoutside')
saveFigureAsEps("prob2_normal_prior.eps", fig);

% Plot Marginal Probabilities
fig = figure("Name", "Marginal Probabilities with Normal Prior");
tl = tiledlayout(2, 1, "Parent", fig);
title(tl, "Marginal Probabilities with Normal Prior")
ax = nexttile(1);
histogram(mout(1,:), "Normalization", "pdf")
title("p_\epsilon(\epsilon)")
xlabel("\epsilon")
ylabel("PDF")
xlim([ep(1) ep(end)])
grid on
grid minor
ax = nexttile(2);
histogram(mout(2,:), "Normalization", "pdf")
title("p_T(T)")
xlabel("T")
ylabel("PDF")
xlim([T(1) T(end)])
grid on
grid minor
saveFigureAsEps("prob2_normal_prior_marginal_prob.eps", fig);

%% Timer

fprintf("Elapsed Time: %6.2f minutes\n\n", toc/60)

