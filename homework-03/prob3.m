close all
clear
clc

addpath(genpath(fullfile("..", "PEIP-master", "Lib")))

% Save figures as *.eps
saveFigureAsEps = @(name, fig)(exportgraphics(fig, fullfile(pwd, "latex", "images", name)));

% Load Provided Data
load(fullfile(pwd, "data", "crosswell.mat"))

%% Problem 3 - Exercise 3 in Section 4.10

fprintf("************************\n")
fprintf("Problem 3\n")
fprintf("************************\n\n")

% Givens
depths = 50 : 50 : 1550;
sigma = 0.5e-3;
GRID_SIZE = 16;

% Size of G
[m, n] = size(G);
fprintf("Model Operator G: m = %d, n = %d\n\n", m, n)

% Visualize G
fig = figure("Name", "Model Operator G");
ax = gca;
hold(ax, "on")
colormap('gray')
imagesc(G)
clim([min(min(G)) max(max(G))])
ax.XTick = (1 : GRID_SIZE : n+1) - 1;
ax.XTick = (1 : GRID_SIZE : n+1) - 1;
ax.YDir = "reverse";
xlabel("columns")
ylabel("rows")
ax.XLim = [ax.XTick(1) - 0.5, ax.XTick(end) + 0.5];
ax.YLim = [ax.YTick(1) - 0.5, ax.YTick(end) + 0.5];
title("Model Operator G")
colorbar(ax, "eastoutside")
saveFigureAsEps("prob3_model_operator_G.eps", fig)

% Set up for weighted least squares
W = sigma * eye(m);
Gw = W * G;
dw = W * dn;


%% Part A - Solve the Inverse Problem using TSVD and L-Curve

fprintf("************************\n")
fprintf("Part A - Using TSVD\n")
fprintf("************************\n\n")

% Rank of G
p = rank(Gw);
fprintf("p = rank(Gw) = %d\n\n", p)

% SVD
[U, S, V] = svd(Gw);

% Create L-Curve for TVSD
modelResiduals = zeros(p, 1);
modelLengths = zeros(p, 1);
for k = 1 : p

    % Truncate
    Up = U(:,1:k);
    Sp = S(1:k,1:k);
    Vp = V(:,1:k);

    % Solve for Model
    m_test = Vp * inv(Sp) * Up.' * dw;

    % Determine Model Length and Residuals
    modelResiduals(k) = norm((Gw*m_test - dw), 2);
    modelLengths(k) = norm(m_test, 2);

end

% Determine Truncation Value
rho = vecnorm([modelResiduals, modelLengths], 2, 2);
[~, pp] = min(rho);

% TEST: Manually pick "pp"
pp = 100;

% Report Truncation Value
fprintf("Selected truncation value: %d\n\n", pp)

% Plot L-Curve
fig = figure("Name", "Part A - TSVD L-Curve");
ax = gca;
hold(ax, "on")
plot(ax, modelResiduals, modelLengths, 'b.')
plot(ax, modelResiduals(pp), modelLengths(pp), 'ro')
ax.XAxis.Scale = "log";
ax.YAxis.Scale = "log";
title("TSVD L-Curve")
xlabel("||Gm - d||_2")
ylabel("||m||_2")
grid on
grid minor
axis equal
saveFigureAsEps("prob3_partA_TSVD_L_Curve.eps", fig)

% Solve the Inverse Problem with Truncated Value
Upp = U(:,1:pp);
Spp = S(1:pp,1:pp);
Vpp = V(:,1:pp);
m_tsvd = Vpp * inv(Spp) * Upp.' * dw;

% Visualize m_tsvd
fig = figure("Name", "Part A - Model via TSVD");
ax = gca;
hold(ax, "on")
colormap('gray')
imagesc(reshape(m_tsvd, [GRID_SIZE GRID_SIZE]).')
clim([min(min(m_tsvd)) max(max(m_tsvd))])
ax.XTick = (1 : GRID_SIZE + 1) - 1;
ax.XTick = (1 : GRID_SIZE + 1) - 1;
ax.YDir = "reverse";
xlabel("columns")
ylabel("rows")
ax.XLim = [ax.XTick(1) - 0.5, ax.XTick(end) + 0.5];
ax.YLim = [ax.YTick(1) - 0.5, ax.YTick(end) + 0.5];
title("m_T_S_V_D")
colorbar(ax, "eastoutside")
saveFigureAsEps("prob3_partA_m_tsvd.eps", fig)

