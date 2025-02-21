%% Problem 1
%
% Exercise 1 from section 2.7

close all
clear
clc

warning off

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
M = length(x);

% Build up the operator G
G = [ones(M, 1), x];

% Solve the least squares problem
m_ls = (G.' * G) \ G.' * t;
t_0 = m_ls(1);
s_2 = m_ls(2);

% Print results
fprintf("Least-Squares Solution\n")
fprintf("\tt_0: %8.6f \n", t_0)
fprintf("\ts_2: %8.6f \n", s_2)

% Compute residuals
r = t - G*m_ls;

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
plot(x, r, 'rx')
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
A = (G.' * G) \ G.';
C = A * A.';

% Print results
fprintf("Model Covariance Matrix\n")
disp(C)

% Comment on the results
% Take note that matrix is positive semi-definite as expected. The model
% parameter t_0 has greater uncertianty than s_2. These two model
% parameters have a slight negative correlation.

% Space out print outs
fprintf("\n\n")


%% Part C

% Print preamble
fprintf("**************************\n")
fprintf("Part C\n")
fprintf("**************************\n\n")



% Space out print outs
fprintf("\n\n")

