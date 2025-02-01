close all
clear
clc

%% Problem 01

% Discritize the integral into 20 equally spaced points
n = 20;
dz = 1/n;
z = 0 : dz : 1;

% Select values of x
nx = 30;  % just picking something different than 20 for now
dx = 1/nx;
x = 0 : dx : 1;

% Create Function Handles
kernal = @(x, z)(5.*sin(x.*z));
m_true = @(x)(10.*x.*sin(x));
fredholm_inside = @(x, z)(kernal(x, z) * m_true(x));
fredhom_ifk = @(x, z)(integral(fredholm_inside(x, z), 0, 1));




% % Plot Results
% fig = figure("Name", "IFK vs. True Model");
% hold on
% plot(x, trueModel(x), 'k', "LineWidth", 2)
% plot(x, d_true(x), 'r.', 'MarkerSize', 5)
% title("IFK vs. True Model")
% xlabel("x")
% ylabel("m_t_r_u_e(x)")
% grid on
% grid minor
% hold off