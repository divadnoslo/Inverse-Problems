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
p = rank(G);
fprintf("p = rank(G) = %d\n\n", p)

% SVD
[U, S, V] = svd(G);
s = diag(S);

% Curious why this doesn't work...
%
% % Create L-Curve for TVSD
% modelResiduals = zeros(p, 1);
% modelLengths = zeros(p, 1);
% for k = 1 : p
% 
%     % Truncate
%     Up = U(:,1:k);
%     Sp = S(1:k,1:k);
%     Vp = V(:,1:k);
% 
%     % Solve for Model
%     m_test = Vp * inv(Sp) * Up.' * dw;
% 
%     % Determine Model Length and Residuals
%     modelResiduals(k) = norm((Gw*m_test - dw), 2);
%     modelLengths(k) = norm(m_test, 2);
% 
% end
% 
% % Determine Truncation Value
% rho = vecnorm([modelResiduals, modelLengths], 2, 2);
% [~, pp] = min(rho);

% Use provided function 'l_curve_tsvd()'
[rho, eta, reg_param] = l_curve_tsvd(U, s, dn);
[~, ireg_corner] = l_curve_corner(rho, eta, reg_param);
rho_corner = rho(ireg_corner);
eta_corner = eta(ireg_corner);

% Determine Truncation Value
pp = ireg_corner;
fprintf("Selected truncation value: %d\n\n", pp)

% Plot L-Curve
fig = figure("Name", "Part A - TSVD L-Curve");
ax = gca;
hold(ax, "on")
plot(ax, rho, eta, 'b.')
plot(ax, rho_corner, eta_corner, 'ro')
ax.XAxis.Scale = "log";
ax.YAxis.Scale = "log";
title("TSVD L-Curve")
xlabel("||Gm - d||_2")
ylabel("||m||_2")
grid on
grid minor
saveFigureAsEps("prob3_partA_TSVD_L_Curve.eps", fig)

% SVD Weight
[U, S, V] = svd(Gw);
s = diag(S);

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


%% Part B - Solve the Inverse Problem using Zeroth-Order Tikhonov Regularization

fprintf("************************\n")
fprintf("Part B - Using Zeroth-Order Tikhonov Regularization\n")
fprintf("************************\n\n")

% SVD
[U, S, V] = svd(G);
s = diag(S);

% L-Curve
[rho, eta, reg_param] = l_curve_tikh_svd(U, s, dn, 1000);
[alpha_tikh, ireg_corner] = l_curve_corner(rho, eta, reg_param);
rho_corner = rho(ireg_corner);
eta_corner = eta(ireg_corner);

% Report alpha_tikh
fprintf("alpha_tikh: %4.3e\n\n", alpha_tikh)

% Plot L-Curve
fig = figure("Name", "Part B - Zeroth-Order Tikhonov L-Curve");
ax = gca;
hold(ax, "on")
plot(ax, rho, eta, 'b', 'LineWidth', 3)
plot(ax, rho_corner, eta_corner, 'rx', 'MarkerSize', 20)
ax.XAxis.Scale = "log";
ax.YAxis.Scale = "log";
title("Zeroth-Order Tikhonov L-Curve")
xlabel("||Gm - d||_2")
ylabel("||m||_2")
grid on
grid minor
saveFigureAsEps("prob3_partB_tikh_L_Curve.eps", fig)

% Compute Zeroth-Order Tikhonov Solution
m_tikh = (Gw.'*Gw + (alpha_tikh^2)*eye(n)) \ Gw.' * dw;

% Visualize m_tikh
fig = figure("Name", "Part B - Model via Zeroth-Order Tikhonov Regularization");
ax = gca;
hold(ax, "on")
colormap('gray')
imagesc(reshape(m_tikh, [GRID_SIZE GRID_SIZE]).')
clim([min(min(m_tikh)) max(max(m_tikh))])
ax.XTick = (1 : GRID_SIZE + 1) - 1;
ax.XTick = (1 : GRID_SIZE + 1) - 1;
ax.YDir = "reverse";
xlabel("columns")
ylabel("rows")
ax.XLim = [ax.XTick(1) - 0.5, ax.XTick(end) + 0.5];
ax.YLim = [ax.YTick(1) - 0.5, ax.YTick(end) + 0.5];
title("Model via Zeroth-Order Tikhonov Regularization")
colorbar(ax, "eastoutside")
saveFigureAsEps("prob3_partB_m_tikh.eps", fig)


%% Part C - Solve the Inverse Problem using Second-Order Tikhonov Regularization

fprintf("************************\n")
fprintf("Part C - Using Second-Order Tikhonov Regularization\n")
fprintf("************************\n\n")

% Calling the "l_code." script from a local function below.
L2 = getL();

% Vizualize L because I am curious
fig = figure("Name", "Roughing Matrix L");
ax = gca;
hold(ax, "on")
colormap('gray')
imagesc(L2)
clim([min(min(L2)) max(max(L2))])
ax.YDir = "reverse";
xlabel("columns")
ylabel("rows")
title("Roughing Matrix L_2")
colorbar(ax, "eastoutside")
saveFigureAsEps("prob3_partC_roughing_matrix_L.eps", fig)

% We didn't cover the GSVD, but it looks like I have to call it?
[U2, V2, X2, LAM2, MU2] = gsvd(G, L2);

% L-Curve
[rho, eta, reg_param, ms] = l_curve_tikh_gsvd(U2, dw, X2, LAM2, MU2, G, L2, 1200);
[alpha_tikh_L2, ireg_corner, kappa] = l_curve_corner(rho, eta, reg_param);
rho_corner = rho(ireg_corner);
eta_corner = eta(ireg_corner);

% Report alpha_tikh
fprintf("alpha_tikh_L2: %4.3e\n\n", alpha_tikh_L2)

% Plot L-Curve
fig = figure("Name", "Part C - Second-Order Tikhonov L-Curve");
ax = gca;
hold(ax, "on")
plot(ax, rho, eta, 'b', 'LineWidth', 3)
plot(ax, rho_corner, eta_corner, 'rx', 'MarkerSize', 20)
ax.XAxis.Scale = "log";
ax.YAxis.Scale = "log";
title("Second-Order Tikhonov L-Curve")
xlabel("||Gm - d||_2")
ylabel("||m||_2")
grid on
grid minor
saveFigureAsEps("prob3_partC_tikh_L_Curve.eps", fig)

% Compute Second-Order Tikhonov Solution
m_tikh_L2 = (Gw.'*Gw + (alpha_tikh_L2^2).*(L2.')*L2) \ Gw.' * dw;

% Visualize m_tikh
fig = figure("Name", "Part C - Model via Second-Order Tikhonov Regularization");
ax = gca;
hold(ax, "on")
colormap('gray')
imagesc(reshape(m_tikh_L2, [GRID_SIZE GRID_SIZE]).')
clim([min(min(m_tikh_L2)) max(max(m_tikh_L2))])
ax.XTick = (1 : GRID_SIZE + 1) - 1;
ax.XTick = (1 : GRID_SIZE + 1) - 1;
ax.YDir = "reverse";
xlabel("columns")
ylabel("rows")
ax.XLim = [ax.XTick(1) - 0.5, ax.XTick(end) + 0.5];
ax.YLim = [ax.YTick(1) - 0.5, ax.YTick(end) + 0.5];
title("Model via Second-Order Tikhonov Regularization")
colorbar(ax, "eastoutside")
saveFigureAsEps("prob3_partC_m_tikh.eps", fig)


%% Local Functions

% Just copying the "l_code.m" into a function to perserve function scope
function L = getL()

    L=zeros(14*14,256);
    k=1;
    for i=2:15
        for j=2:15
            M=zeros(16,16);
            M(i,j)=-4;
            M(i,j+1)=1;
            M(i,j-1)=1;
            M(i+1,j)=1;
            M(i-1,j)=1;
            L(k,:)=reshape(M,256,1)';
            k=k+1;
        end
    end

end