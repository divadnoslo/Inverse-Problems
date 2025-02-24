%% Problem 2
%
% Exercise 5 from section 2.7

close all
clear
clc

warning off

% Add "Lib" to serach path
addpath(genpath(fullfile("..", "PEIP-master", "Lib")))

% Save figures as *.eps
saveFigureAsEps = @(name, fig)(exportgraphics(fig, fullfile(pwd, "latex", "images", name)));

% Print Preamble
fprintf("Homework 2, Problem 2\n\n")


%% Homework Solution

% Set up sample points
x = (-1 : 0.1 : 1).';
y = (-1 : 0.1 : 1).';

% Build up operator G
m = length(x);
n = 20;
G = zeros(m, n);
for ii = 1 : m
    for jj = 1 : n
        G(ii,jj) = x(ii)^(jj - 1);
    end
end

% Compute least squares solution
m_L2 = (G.' * G) \ G.' * y;

% Display model parameters
disp("L2 Regression Model Parameters")
format shorte
disp(m_L2)
format shortg

% Compute data for comparison
xs = (-1 : 1e-3 : 1).';
y_true = xs;
y_fit = zeros(length(xs), 1);
for k = 1 : length(xs)
    for kk = 1 : n
        y_fit(k) = y_fit(k) + m_L2(kk) * xs(k)^(kk - 1);
    end
end

% Plot comparison
fig = figure("Name", "Truth vs. Least-Squares Fit");
tl = tiledlayout(2, 1, "Parent", fig);
title(tl, "Homework 2, Problem 2")

ax = nexttile(1);
hold(ax, "on")
plot(xs, y_true, 'k')
plot(xs, y_fit, 'b')
plot(x, y, 'g.', 'MarkerSize', 10)
title("Truth vs. L2 Fit")
xlabel("x")
ylabel("y")
grid on
grid minor
legend(["Truth", "L2 Fit", "Data Points"], "Location", "eastoutside")

ax = nexttile(2);
hold(ax, "on")
plot(xs, y_fit - y_true, 'r')
plot(x, y - G*m_L2, 'm.', 'MarkerSize', 10)
title("Residuals")
xlabel("x")
ylabel("\Delta y")
grid on
grid minor
legend(["Fit - Truth", "L2 Fit Residuals"], "Location", "eastoutside")

linkaxes(tl.Children, 'x')
saveFigureAsEps("prob2_1.eps", fig)