close all
clear
clc

% Add to Path
addpath(genpath(fullfile("..", "PEIP-master", "Lib")))

% Save figures as *.eps
saveFigureAsEps = @(name, fig)(exportgraphics(fig, fullfile(pwd, "latex", "images", name)));

% Load Data
load(fullfile(pwd, "data", "be10.mat"))


%% Problem 2 - Exercise 8 in Section 9.6

% Givens
lambda = 4.998e-7;
P0 = 3.2;
mu = 0.0166;
z = depths;

% Function Handle
P = @(z)(P0 * exp(-mu*z));


%% Part A

% Points
ep = 5e-6 : 1e-5 : 1e-3;
T = 500 : 1000 : 199500; 

% Make Grid
[X, Y] = meshgrid(ep, T);

% Compute Chi-Squared over the Region
Z = zeros(size(X));
[m, n] = size(X);
for ii = 1 : m
    for jj = 1 : n
        Z(ii,jj) = norm(prob2Function([X(ii,jj); Y(ii,jj)]))^2;
    end
end  

% Find the minimum
[minValues, rowIndices] = min(Z);
[globalMin, colIndex] = min(minValues);
rowIndex = rowIndices(colIndex);

% Plot the resulting surface
fig = figure("Name", "Chi-Squared Surface");
ax = gca;
hold(ax, "on")
sh = surf(X, Y, Z);
ph = plot3(ep(rowIndex), T(colIndex), globalMin, 'r.', 'MarkerSize', 25);
title("\chi^2 Surface")
xlabel("\epsilon")
ylabel("T")
zlabel("\chi^2")
grid on
view(75, 35)
legend(ph, "Minimum", "Location", "eastoutside")
saveFigureAsEps("prob2_chi_2_surface.eps", fig);

% Print Results
fprintf("Finding the Minimum of the Chi-Squared Surface\n")
fprintf("\tValue of Epsilon: %6.6e\n", ep(rowIndex))
fprintf("\tValue of T:       %6.6e\n", T(rowIndex))
fprintf("\tMinimum Value for Chi-Squared: %6.3f\n", globalMin)


%% Part B

