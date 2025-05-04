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
ep = 5e-6 : 1e-5 : 1e-3;
T = 500 : 1000 : 199500; 

% Make Grid
[X, Y] = meshgrid(ep, T);

% Compute Chi-Squared over the Region
chi2 = zeros(size(X));
[m, n] = size(X);
for ii = 1 : m
    for jj = 1 : n
        chi2(ii,jj) = norm(model2residuals([X(ii,jj); Y(ii,jj)]))^2;
    end
end  


%% Part A

