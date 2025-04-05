close all
clear
clc

% Add to Path
addpath(genpath(fullfile("..", "PEIP-master", "Lib")))
addpath(genpath(fullfile(pwd, "prob2files")))

% Save figures as *.eps
saveFigureAsEps = @(name, fig)(exportgraphics(fig, fullfile(pwd, "latex", "images", name)));


%% Problem 2 - Exercise 6 in Section 7.8

% Start Timer
tic

% Load the raw image
img = double(imread(fullfile(pwd, "prob2files", "prob2files", "inpaint.png")));

% Plot Raw Image
fig = figure("Name", "Original Image");
imagesc(img);
title("Original Image")
colormap(gray);
ax = gca;
set(ax, 'XTick', []);
set(ax, 'YTick', []);
axis equal
saveFigureAsEps("prob2_original_image.eps", fig)

% Convert image to vector
m = reshape(img, 512^2, 1);

% Mask Good and Bad Pixels
maskBad = m == 255;
maskGood = ~maskBad;

% Get Data
d = m;
d(maskBad) = 0;



% Stop Timer
fprintf("Elapsed Minutes: %4.2f minutes\n\n", toc/60)