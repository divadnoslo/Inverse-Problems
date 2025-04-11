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
[p, ~] = size(img);
m = reshape(img, p^2, 1);

% Mask Good and Bad Pixels
maskBad = m == 255;
maskGood = ~maskBad;

% Get Model Operator
I = speye(p^2);
I_good = I;
I_good(maskBad,:) = [];

% Get Data
d = m;
d(maskBad) = [];

% Create L Matrix
L = -2 * speye(p^2, p^2 + p);
I_1 = [sparse(p^2,1), speye(p^2, p^2 + p - 1)];
I_2 = [sparse(p^2,p), speye(p^2)];
L = L + I_1 + I_2;
L = L(1:p^2, 1:p^2);
rowsToRemove = mod(1:p^2, p) == 0;
L(rowsToRemove,:) = [];

% Regularize
alpha = 1;
[m_new, bestObj, iter, objectiveValues] = admml1reg(I_good, d, L, alpha);

% Plot Result
fig = figure("Name", "New Image");
imagesc(reshape(m_new, p, p));
title(sprintf("Regularized Image: alpha = %4.2e", alpha))
colormap(gray);
ax = gca;
set(ax, 'XTick', []);
set(ax, 'YTick', []);
axis equal
saveFigureAsEps("prob2_regularized_image.eps", fig)

% Regularize with Low Alpha
alpha = 0.1;
[m_new, bestObj, iter, objectiveValues] = admml1reg(I_good, d, L, alpha);

% Plot Result
fig = figure("Name", "New Image - Low Alpha");
imagesc(reshape(m_new, p, p));
title(sprintf("Regularized Image: alpha = %4.2e", alpha))
colormap(gray);
ax = gca;
set(ax, 'XTick', []);
set(ax, 'YTick', []);
axis equal
saveFigureAsEps("prob2_regularized_image_low_alpha.eps", fig)

% Regularize with High Alpha
alpha = 10;
[m_new, bestObj, iter, objectiveValues] = admml1reg(I_good, d, L, alpha);

% Plot Result
fig = figure("Name", "New Image - High Alpha");
imagesc(reshape(m_new, p, p));
title(sprintf("Regularized Image: alpha = %4.2e", alpha))
colormap(gray);
ax = gca;
set(ax, 'XTick', []);
set(ax, 'YTick', []);
axis equal
saveFigureAsEps("prob2_regularized_image_high_alpha.eps", fig)

% Stop Timer
fprintf("Elapsed Minutes: %4.2f minutes\n\n", toc/60)

