close all
clear
clc

% Add to Path
addpath(genpath(fullfile("..", "PEIP-master", "Lib")))
addpath(genpath(fullfile(pwd, "prob1files")))

% Save figures as *.eps
saveFigureAsEps = @(name, fig)(exportgraphics(fig, fullfile(pwd, "latex", "images", name)));


%% Problem 1 - Modifying Code from Example 6.2

% Start Timer
tic

% Seed RNG
rand('state', 0);
randn('state', 0);

% Load the raw image
img = double(imread(fullfile(pwd, "prob1files", "image.png")));

% Build the G matrix
G = blur(200, 3, 15);
[m, n] = size(G);

% Compute the blurred image
d = G*reshape(img, 40000, 1);

% Add noise as appropriate
dn = d + 2.0e-2*randn(size(d));

% Plot Raw Image
fig = figure("Name", "Raw Image");
imagesc(img);
title("Raw Image")
colormap(gray);
ax = gca;
set(ax, 'XTick', []);
set(ax, 'YTick', []);
axis equal
saveFigureAsEps("prob1_raw_image.eps", fig)

% Plot Noisy Image
fig = figure("Name", "Noisy Image");
imagesc(reshape(dn, 200, 200));
title("Noisy Image")
colormap(gray);
ax = gca;
set(ax, 'XTick', []);
set(ax, 'YTick', []);
axis equal
saveFigureAsEps("prob1_noisy_image.eps", fig)

% Use Landweber Iterations for Method of Steepest Decent
m0 = zeros(n, 1);
omega = 0.95 * (2/(svds(G,1)^2));
k = 500;
X = landweber(G, m0, dn, omega, k);

%% Evaluate Iterations

e = [10, 50, 100, 150, 250, 500];
assert(e(end) <= k)
for iter = e
    fig = figure("Name", sprintf("%d Iterations Deblurred Image", iter));
    imagesc(reshape(X(:,iter), 200, 200));
    title(sprintf("Deblurred Image - %d Iterations", iter))
    colormap(gray);
    ax = gca;
    set(ax, 'XTick', []);
    set(ax, 'YTick', []);
    axis equal
    saveFigureAsEps(sprintf("prob1_deblurred_image_%d_iterations.eps", iter), fig)
end

% Stop Timer
fprintf("Elapsed Minutes: %4.2f minutes\n\n", toc/60)

