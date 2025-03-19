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
fprintf("---B---\n")
[U, S, V] = svd(G);
Up = U(:,1:p);
Sp = S(1:p,1:p);
Vp = V(:,1:p);
U0 = U(:,p+1:end);
V0 = V(:,p+1:end);

fprintf("Size of Up: ")
disp(size(Up))
fprintf("Size of Vp: ")
disp(size(Vp))
fprintf("Size of U0: ")
disp(size(U0))
fprintf("Size of V0: ")
disp(size(V0))

fig = figure("Name", "Null Space Examples");
tl = tiledlayout(1, 2, "Parent", fig);

ax = nexttile(1);
hold(ax, "on")
stem(1:m, U0)
title("Data Null Space - U_0 = U_._,_3_2")
xlabel("Element")
ylabel("Element Value")
grid on
grid minor

ax = nexttile(2);
hold(ax, "on")
colormap('gray')
imagesc(reshape(V0(:,1), [GRID_SIZE GRID_SIZE]).')
clim([min(V0(:,1)) max(V0(:,1))])
ax.XTick = 0 : GRID_SIZE;
ax.YTick = 0 : GRID_SIZE;
ax.YDir = "reverse";
xlabel("columns")
ylabel("rows")
ax.XLim = [ax.XTick(1) - 0.5, ax.XTick(end) + 0.5];
ax.YLim = [ax.YTick(1) - 0.5, ax.YTick(end) + 0.5];
title("Model Null Space - V_._,_3_2")
axis equal
colorbar(ax, "eastoutside")

saveFigureAsEps("prob2_partA_null_space_examples.eps", fig)

% Model Resolution Matrix
Rm = Vp * Vp.';

fig = figure("Name", "Model Resolution Matrix");
tl = tiledlayout(1, 2, "Parent", fig);

ax = nexttile(1);
hold(ax, "on")
colormap('gray')
imagesc(Rm)
clim([-0.15 0.15])
ax.XTick = 1 : GRID_SIZE : GRID_SIZE^2;
ax.YTick = 1 : GRID_SIZE : GRID_SIZE^2;
ax.YDir = "reverse";
xlabel("columns")
ylabel("rows")
ax.XLim = [ax.XTick(1) - 0.5, ax.XTick(end) + 0.5];
ax.YLim = [ax.YTick(1) - 0.5, ax.YTick(end) + 0.5];
title("Model Resolution Matrix")
axis equal
colorbar(ax, "eastoutside")

Rm_diag = reshape(diag(Rm), [GRID_SIZE GRID_SIZE]).';

ax = nexttile(2);
hold(ax, "on")
colormap('gray')
imagesc(Rm_diag)
clim([-0.15 0.15])
ax.XTick = 0 : GRID_SIZE;
ax.YTick = 0 : GRID_SIZE;
ax.YDir = "reverse";
xlabel("columns")
ylabel("rows")
ax.XLim = [ax.XTick(1) - 0.5, ax.XTick(end) + 0.5];
ax.YLim = [ax.YTick(1) - 0.5, ax.YTick(end) + 0.5];
title("Model Resolution Matrix - Diagonal Elements")
axis equal
colorbar(ax, "eastoutside")

saveFigureAsEps("prob2_partA_model_resolution_matrix.eps", fig)

% Part C - Compute the Model
fprintf("---C---\n")
m_dagger = pinv(G) * d;

[maxSlowness, maxSlownessInd] = max(m_dagger);
[minSlowness, minSlownessInd] = min(m_dagger);

maxSlownessY = floor(maxSlownessInd / GRID_SIZE) + 1;
maxSlownessX = mod(maxSlownessInd, GRID_SIZE);

minSlownessY = floor(minSlownessInd / GRID_SIZE) + 1;
minSlownessX = mod(minSlownessInd, GRID_SIZE);

fig = figure("Name", "Estimated Model Parameters");
ax = gca;
hold(ax, "on")
colormap('gray')
imagesc(reshape(m_dagger, [GRID_SIZE GRID_SIZE]).')
plot(maxSlownessX, maxSlownessY, 'rx', 'MarkerSize', 7)
plot(minSlownessX, minSlownessY, 'ro', 'MarkerSize', 7)
clim([min(m_dagger) max(m_dagger)])
ax.XTick = 0 : GRID_SIZE;
ax.YTick = 0 : GRID_SIZE;
ax.YDir = "reverse";
xlabel("columns")
ylabel("rows")
ax.XLim = [ax.XTick(1) - 0.5, ax.XTick(end) + 0.5];
ax.YLim = [ax.YTick(1) - 0.5, ax.YTick(end) + 0.5];
title("Estimated Model Parameters")
axis equal
colorbar(ax, "eastoutside")
legend(["Max Slowness", "Min Slowness"], "Location", "westoutside")

saveFigureAsEps("prob2_partA_estimated_model_parameters.eps", fig)

% Print Maximums and Minimums
fprintf("Maximum Slowness Value: %8.3e sec/m\n", maxSlowness)
fprintf("Minimum Slowwnes Value: %8.3e sec/m\n", minSlowness)
fprintf("\n")

% Interpret Siesmic Velocity
SV = VELOCITY + (1./m_dagger);

[maxVelocity, maxVelocityInd] = max(SV);
[minVelocity, minVelocityInd] = min(SV);

maxVelocityY = floor(maxVelocityInd / GRID_SIZE) + 1;
maxVelocityX = mod(maxVelocityInd, GRID_SIZE);

minVelocityY = floor(minVelocityInd / GRID_SIZE) + 1;
minVelocityX = mod(minVelocityInd, GRID_SIZE);

fig = figure("Name", "Interpreted Seismic Velocity");
ax = gca;
hold(ax, "on")
colormap('gray')
imagesc(reshape(SV, [GRID_SIZE GRID_SIZE]).')
plot(maxVelocityX, maxVelocityY, 'rx', 'MarkerSize', 7)
plot(minVelocityX, minVelocityY, 'ro', 'MarkerSize', 7)
clim([min(SV) max(SV)])
ax.XTick = 0 : GRID_SIZE;
ax.YTick = 0 : GRID_SIZE;
ax.YDir = "reverse";
xlabel("columns")
ylabel("rows")
ax.XLim = [ax.XTick(1) - 0.5, ax.XTick(end) + 0.5];
ax.YLim = [ax.YTick(1) - 0.5, ax.YTick(end) + 0.5];
title("Interpreted Siesmic Velocity")
axis equal
colorbar(ax, "eastoutside")
legend(["Max Velocity", "Min Velocity"], "Location", "westoutside")

saveFigureAsEps("prob2_partA_seismic_velocity.eps", fig)

% Print Maximums and Minimums
fprintf("Maximum Velocity Value: %8.3f m/sec\n", maxVelocity)
fprintf("Minimum Velocity Value: %8.3f m/sec\n", minVelocity)
fprintf("\n")

% Part D - Wild Model
fprintf("---D---\n")

m_wild = 0.25*V(:,32) + 0.50*V(:,127) + 0.25*V(:,183);
d_wild = G * m_wild;

fig = figure("Name", "Wild Model and Resulting Data");
tl = tiledlayout(2, 1, "Parent", fig);

ax = nexttile(1);
hold(ax, "on")
stem(1:n, m_wild)
title("m_w_i_l_d")
xlabel("Model Vector Index")
ylabel("Model Vector Values")
grid on
grid minor

ax = nexttile(2);
hold(ax, "on")
stem(1:m, d_wild)
title("d_w_i_l_d")
xlabel("Data Vector Index")
ylabel("Data Vector Values")
grid on
grid minor

saveFigureAsEps("prob2_partA_wild_model_and_data.eps", fig)


%% Adding in Diagonal Scanes

disp("******************************")
disp("Part B")
disp("******************************")
fprintf("\n")

% Load Data (Format: x1, y1, x2, y2, t)
load(fullfile(pwd, "data", "diag1scan.mat"));
load(fullfile(pwd, "data", "diag2scan.mat"));

[numDiag1Scans, ~] = size(diag1scan);
[numDiag2Scans, ~] = size(diag2scan);

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
for k = 1 : numDiag1Scans
    quiver(...
        ax, ...
        diag1scan(k,1), diag1scan(k,2), ...
        diag1scan(k,3) - diag1scan(k,1), diag1scan(k,4) - diag1scan(k,2), ...
        'Color', 'g', ...
        'AutoScaleFactor', 1, ...
        'MaxHeadSize', 0.05)
end
for k = 1 : numDiag2Scans
    quiver(...
        ax, ...
        diag2scan(k,1), diag2scan(k,2), ...
        diag2scan(k,3) - diag2scan(k,1), diag2scan(k,4) - diag2scan(k,2), ...
        'Color', 'k', ...
        'AutoScaleFactor', 1, ...
        'MaxHeadSize', 0.05)
end
title("Adding in Diagonal Scans")
xlabel("East/West [m]")
ylabel("North/South [m]")
saveFigureAsEps("prob2_partB_scans_vizualization.eps", fig)

% Create d and G
m = 32 + numDiag1Scans + numDiag2Scans;
n = 256;
d = [rowscan(:,5); colscan(:,5); diag1scan(:,5); diag2scan(:,5)];
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

% Each Diagonal 1 Scan
for ii = 1 : GRID_SIZE
    G(32 + (ii-1) + (1:GRID_SIZE), (ii-1)*GRID_SIZE + (1:GRID_SIZE)) = sqrt(2) * eye(GRID_SIZE);
end

% Each Diagonal 2 Scan
Ir = flip(eye(GRID_SIZE), 1);
for ii = 1 : GRID_SIZE
    G(32 + 31 + (ii-1) + (1:GRID_SIZE), (ii-1)*GRID_SIZE + (1:GRID_SIZE)) = sqrt(2) * Ir;
end

% Vizualize and Confirm G
fig = figure("Name", "Model Operator G");
ax = gca;
hold(ax, "on")
colormap('gray')
imagesc(G)
clim([0 sqrt(2)])
ax.XTick = (1 : GRID_SIZE : n+1) - 1;
ax.YTick = 1 : m;
ax.YDir = "reverse";
xlabel("columns")
ylabel("rows")
ax.XLim = [ax.XTick(1) - 0.5, ax.XTick(end) + 0.5];
ax.YLim = [ax.YTick(1) - 0.5, ax.YTick(end) + 0.5];
title("Model Operator G")
colorbar(ax, "eastoutside")
axis equal
saveFigureAsEps("prob2_partB_model_operator_G.eps", fig)

% Part A - Rank of G
fprintf("---A---\n")
p = rank(G);
fprintf("rank(G) = %d\n\n", p)

% Part B - Lots of Stuff
fprintf("---B---\n")
[U, S, V] = svd(G);
Up = U(:,1:p);
Sp = S(1:p,1:p);
Vp = V(:,1:p);
U0 = U(:,p+1:end);
V0 = V(:,p+1:end);

fprintf("Size of Up: ")
disp(size(Up))
fprintf("Size of Vp: ")
disp(size(Vp))
fprintf("Size of U0: ")
disp(size(U0))
fprintf("Size of V0: ")
disp(size(V0))

fig = figure("Name", "Null Space Examples");
tl = tiledlayout(1, 2, "Parent", fig);

ax = nexttile(1);
hold(ax, "on")
stem(1:m, U0(:,1))
title("Data Null Space - U_0 = U_._,_8_8")
xlabel("Element")
ylabel("Element Value")
grid on
grid minor

ax = nexttile(2);
hold(ax, "on")
colormap('gray')
imagesc(reshape(V0(:,1), [GRID_SIZE GRID_SIZE]).')
clim([min(V0(:,1)) max(V0(:,1))])
ax.XTick = 0 : GRID_SIZE;
ax.YTick = 0 : GRID_SIZE;
ax.YDir = "reverse";
xlabel("columns")
ylabel("rows")
ax.XLim = [ax.XTick(1) - 0.5, ax.XTick(end) + 0.5];
ax.YLim = [ax.YTick(1) - 0.5, ax.YTick(end) + 0.5];
title("Model Null Space - V_._,_8_8")
axis equal
colorbar(ax, "eastoutside")

saveFigureAsEps("prob2_partB_null_space_examples.eps", fig)

% Model Resolution Matrix
Rm = Vp * Vp.';

fig = figure("Name", "Model Resolution Matrix");
tl = tiledlayout(1, 2, "Parent", fig);

ax = nexttile(1);
hold(ax, "on")
colormap('gray')
imagesc(Rm)
clim([min(min(Rm)) max(max(Rm))])
ax.XTick = 1 : GRID_SIZE : GRID_SIZE^2;
ax.YTick = 1 : GRID_SIZE : GRID_SIZE^2;
ax.YDir = "reverse";
xlabel("columns")
ylabel("rows")
ax.XLim = [ax.XTick(1) - 0.5, ax.XTick(end) + 0.5];
ax.YLim = [ax.YTick(1) - 0.5, ax.YTick(end) + 0.5];
title("Model Resolution Matrix")
axis equal
colorbar(ax, "eastoutside")

Rm_diag = reshape(diag(Rm), [GRID_SIZE GRID_SIZE]).';

ax = nexttile(2);
hold(ax, "on")
colormap('gray')
imagesc(Rm_diag)
clim([min(min(Rm)) max(max(Rm))])
ax.XTick = 0 : GRID_SIZE;
ax.YTick = 0 : GRID_SIZE;
ax.YDir = "reverse";
xlabel("columns")
ylabel("rows")
ax.XLim = [ax.XTick(1) - 0.5, ax.XTick(end) + 0.5];
ax.YLim = [ax.YTick(1) - 0.5, ax.YTick(end) + 0.5];
title("Model Resolution Matrix - Diagonal Elements")
axis equal
colorbar(ax, "eastoutside")

saveFigureAsEps("prob2_partB_model_resolution_matrix.eps", fig)

% Part C - Compute the Model
fprintf("---C---\n")
m_dagger = pinv(G) * d;

[maxSlowness, maxSlownessInd] = max(m_dagger);
[minSlowness, minSlownessInd] = min(m_dagger);

maxSlownessY = floor(maxSlownessInd / GRID_SIZE) + 1;
maxSlownessX = mod(maxSlownessInd, GRID_SIZE);

minSlownessY = floor(minSlownessInd / GRID_SIZE) + 1;
minSlownessX = mod(minSlownessInd, GRID_SIZE);

fig = figure("Name", "Estimated Model Parameters");
ax = gca;
hold(ax, "on")
colormap('gray')
imagesc(reshape(m_dagger, [GRID_SIZE GRID_SIZE]).')
plot(maxSlownessX, maxSlownessY, 'rx', 'MarkerSize', 7)
plot(minSlownessX, minSlownessY, 'ro', 'MarkerSize', 7)
clim([min(m_dagger) max(m_dagger)])
ax.XTick = 0 : GRID_SIZE;
ax.YTick = 0 : GRID_SIZE;
ax.YDir = "reverse";
xlabel("columns")
ylabel("rows")
ax.XLim = [ax.XTick(1) - 0.5, ax.XTick(end) + 0.5];
ax.YLim = [ax.YTick(1) - 0.5, ax.YTick(end) + 0.5];
title("Estimated Model Parameters")
axis equal
colorbar(ax, "eastoutside")
legend(["Max Slowness", "Min Slowness"], "Location", "westoutside")

saveFigureAsEps("prob2_partB_estimated_model_parameters.eps", fig)

% Print Maximums and Minimums
fprintf("Maximum Slowness Value: %8.3e sec/m\n", maxSlowness)
fprintf("Minimum Slowwnes Value: %8.3e sec/m\n", minSlowness)
fprintf("\n")

% Interpret Siesmic Velocity
SV = VELOCITY + (1./m_dagger);

[maxVelocity, maxVelocityInd] = max(SV);
[minVelocity, minVelocityInd] = min(SV);

maxVelocityY = floor(maxVelocityInd / GRID_SIZE) + 1;
maxVelocityX = mod(maxVelocityInd, GRID_SIZE);

minVelocityY = floor(minVelocityInd / GRID_SIZE) + 1;
minVelocityX = mod(minVelocityInd, GRID_SIZE);

fig = figure("Name", "Interpreted Seismic Velocity");
ax = gca;
hold(ax, "on")
colormap('gray')
imagesc(reshape(SV, [GRID_SIZE GRID_SIZE]).')
plot(maxVelocityX, maxVelocityY, 'rx', 'MarkerSize', 7)
plot(minVelocityX, minVelocityY, 'ro', 'MarkerSize', 7)
clim([min(SV) max(SV)])
ax.XTick = 0 : GRID_SIZE;
ax.YTick = 0 : GRID_SIZE;
ax.YDir = "reverse";
xlabel("columns")
ylabel("rows")
ax.XLim = [ax.XTick(1) - 0.5, ax.XTick(end) + 0.5];
ax.YLim = [ax.YTick(1) - 0.5, ax.YTick(end) + 0.5];
title("Interpreted Siesmic Velocity")
axis equal
colorbar(ax, "eastoutside")
legend(["Max Velocity", "Min Velocity"], "Location", "westoutside")

saveFigureAsEps("prob2_partB_seismic_velocity.eps", fig)

% Print Maximums and Minimums
fprintf("Maximum Velocity Value: %8.3f m/sec\n", maxVelocity)
fprintf("Minimum Velocity Value: %8.3f m/sec\n", minVelocity)
fprintf("\n")

% Part D - Wild Model
fprintf("---D---\n")

m_wild = 0.25*V(:,88) + 0.50*V(:,132) + 0.25*V(:,201);
d_wild = G * m_wild;

fig = figure("Name", "Wild Model and Resulting Data");
tl = tiledlayout(2, 1, "Parent", fig);

ax = nexttile(1);
hold(ax, "on")
stem(1:n, m_wild)
title("m_w_i_l_d")
xlabel("Model Vector Index")
ylabel("Model Vector Values")
grid on
grid minor

ax = nexttile(2);
hold(ax, "on")
stem(1:m, d_wild)
title("d_w_i_l_d")
xlabel("Data Vector Index")
ylabel("Data Vector Values")
grid on
grid minor

saveFigureAsEps("prob2_partB_wild_model_and_data.eps", fig)

