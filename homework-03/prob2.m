close all
clear
clc

addpath(genpath(fullfile("..", "PEIP-master", "Lib")))

% Save figures as *.eps
saveFigureAsEps = @(name, fig)(exportgraphics(fig, fullfile(pwd, "latex", "images", name)));


%% Problem2

% Problem given constants
VELOCITY = 3000;
SIGMA = 1e-5;
GRID_SIZE = 16;


%% Row and Column Scans Only

disp("******************************")
disp("Part A")
disp("******************************")
fprintf("\n")

% Load Data (Format: x1, y1, x2, y2, t)
load(fullfile(pwd, "data", "rowscan.mat"));
load(fullfile(pwd, "data", "colscan.mat"));

[numRowScans, ~] = size(rowscan);
[numColScans, ~] = size(colscan);

% Vizualize
fig = figure("Name", "Data Vizualization");
ax = gca;
hold(ax, "on")
ax.YDir = "reverse";
ax.XLim = [0 GRID_SIZE];
ax.YLim = [0 GRID_SIZE];
grid(ax, "on")
grid(ax, "minor")
ax.XTick = 0 : GRID_SIZE;
ax.YTick = 0 : GRID_SIZE;
ax.GridAlpha = 1;
for k = 1 : numRowScans
    quiver(...
        ax, ...
        rowscan(k,1), rowscan(k,2), ...
        rowscan(k,3) - rowscan(k,1), rowscan(k,4) - rowscan(k,2), ...
        'Color', 'b', ...
        'AutoScaleFactor', 1, ...
        'MaxHeadSize', 0.05)
end
for k = 1 : numColScans
    quiver(...
        ax, ...
        colscan(k,1), colscan(k,2), ...
        colscan(k,3) - colscan(k,1), colscan(k,4) - colscan(k,2), ...
        'Color', 'r', ...
        'AutoScaleFactor', 1, ...
        'MaxHeadSize', 0.05)
end
title("Row and Column Scans")
xlabel("East/West [m]")
ylabel("North/South [m]")
saveFigureAsEps("prob2_partA_scans_vizualization.eps", fig)

% Create d and G
m = 32;
n = 256;
d = [rowscan(:,5); colscan(:,5)];
assert(length(d) == m)
G = zeros(m, n);

% Each Row Scan
for ii = 1 : GRID_SIZE
        G(ii, (1:GRID_SIZE)+GRID_SIZE*(ii-1)) = ones(GRID_SIZE,1);
end

% Each Column Scan
for jj = 1 : GRID_SIZE
    G(1+GRID_SIZE : 2*GRID_SIZE, (1 + (jj-1)*GRID_SIZE) : jj*GRID_SIZE) = eye(GRID_SIZE);
end

% Vizualize and Confirm G
fig = figure("Name", "Model Operator G");
ax = gca;
hold(ax, "on")
colormap('gray')
imagesc(G)
clim([0 1])
ax.XTick = (1 : GRID_SIZE : n+1) - 1;
ax.YTick = 1 : m;
ax.YDir = "reverse";
xlabel("columns")
ylabel("rows")
ax.XLim = [ax.XTick(1) - 0.5, ax.XTick(end) + 0.5];
ax.YLim = [ax.YTick(1) - 0.5, ax.YTick(end) + 0.5];
title("Model Operator G")
colorbar(ax, "eastoutside")
saveFigureAsEps("prob2_partA_model_operator_G.eps", fig)

% Part A - Rank of G
fprintf("---A---\n")
p = rank(G);
fprintf("rank(G) = %d\n\n", p)

% Part B - Lots of Stuff
[U, S, V] = svd(G);
Up = U(1:p,1:p);
Sp = S(1:p,1:p);
Vp = V(:,1:p);

m_dagger = pinv(G) * d;