close all
clear
clc

addpath(genpath(fullfile("..", "PEIP-master", "Lib")))


%% Let's go through Example 3.1 first...

fprintf("**************************\n")
disp("Example 3.1")
fprintf("**************************\n\n")

% Copy G from the textbook
G = [...
    1 0 0 1 0 0 1 0 0; ...
    0 1 0 0 1 0 0 1 0; ...
    0 0 1 0 0 1 0 0 1; ...
    1 1 1 0 0 0 0 0 0; ...
    0 0 0 1 1 1 0 0 0; ...
    0 0 0 0 0 0 1 1 1; ...
    sqrt(2) 0 0 0 sqrt(2) 0 0 0 sqrt(2); ...
    0 0 0 0 0 0 0 0 sqrt(2)];
[m, n] = size(G);
assert(m == 8);
assert(n == 9);

% SVD
[U, S, V] = svd(G);
disp("Singular Values of G")
s = diag(S);
disp(s)

% Rank of G
p = rank(G);
fprintf("rank(G) = %d\n\n", p)

% Ratio of Singular Values
s_p = s(s > 1e-6);
ratio = s_p(1) / s_p(end);
assert(ratio < p)
fprintf("Ratio vs. Rank: %4.3f < %d (good)\n\n", ratio, p)

% Model Null Space
V0 = V(:, (1:n > p));
disp("Model Null Space (V0)")
disp(V0)

% Visualize
v8 = reshape(V0(:,1), 3, 3).';
v9 = reshape(V0(:,2), 3, 3).';

fig = figure("Name", "Model Null Space");
tl = tiledlayout(1, 2, "Parent", fig);
title(tl, "Model Null Space")

ax = nexttile(1);
hold(ax, "on")
colormap('gray')
imagesc(v8)
clim([-0.5 0.5])
ax.XTick = 1 : 3;
ax.YTick = 1 : 3;
ax.YDir = "reverse";
xlabel("j")
ylabel("i")
ax.XLim = [ax.XTick(1) - 0.5, ax.XTick(end) + 0.5];
ax.YLim = [ax.YTick(1) - 0.5, ax.YTick(end) + 0.5];
axis equal
title("Reshaped Null Space V(:,8)")
colorbar(ax, "eastoutside")

ax = nexttile(2);
hold(ax, "on")
colormap('gray')
imagesc(v9)
clim([-0.5 0.5])
ax.XTick = 1 : 3;
ax.YTick = 1 : 3;
ax.YDir = "reverse";
xlabel("j")
ylabel("i")
ax.XLim = [ax.XTick(1) - 0.5, ax.XTick(end) + 0.5];
ax.YLim = [ax.YTick(1) - 0.5, ax.YTick(end) + 0.5];
axis equal
title("Reshaped Null Space V(:,9)")
colorbar(ax, "eastoutside")

% Data Null Space
U0 = U(:, (1:m > p));
disp("Data Null Space (U0)")
disp(U0)

% Model Resolution Matrix
Vp = V(:,1:p);
Rm = Vp * Vp.';
rmDiag = diag(Rm);
rm = reshape(rmDiag, 3, 3).';

% Vizualize Model Resolution Matrix
fig = figure("Name", "Model Resolution Matrix");
tl = tiledlayout(1, 2, "Parent", fig);
title(tl, "Model Resolution Matrix")

ax = nexttile(1);
hold(ax, "on")
colormap('gray')
imagesc(Rm)
clim([0 1])
ax.XTick = 1 : 9;
ax.YTick = 1 : 9;
ax.YDir = "reverse";
xlabel("j")
ylabel("i")
ax.XLim = [ax.XTick(1) - 0.5, ax.XTick(end) + 0.5];
ax.YLim = [ax.YTick(1) - 0.5, ax.YTick(end) + 0.5];
axis equal
title("Model Resolution Matrix (Rm)")
colorbar(ax, "eastoutside")

ax = nexttile(2);
hold(ax, "on")
colormap('gray')
imagesc(rm)
clim([0 1])
ax.XTick = 1 : 3;
ax.YTick = 1 : 3;
ax.YDir = "reverse";
xlabel("j")
ylabel("i")
ax.XLim = [ax.XTick(1) - 0.5, ax.XTick(end) + 0.5];
ax.YLim = [ax.YTick(1) - 0.5, ax.YTick(end) + 0.5];
axis equal
title("Model Resolution Matrix Diagonals")
colorbar(ax, "eastoutside")

% Set Up Spike Test
m_test = zeros(n, 1);
m_test(5) = 1;
d_test = G * m_test;
disp("Model Spike")
disp(m_test)
disp("Predicted Data")
disp(d_test)

% Conduct Spike Test
m_dagger = pinv(G) * d_test;
disp("Recovered Model from Spike Test")
disp(m_dagger)

% Vizualize
fig = figure("Name", "Spike Test");
tl = tiledlayout(1, 2, "Parent", fig);
title(tl, "Spike Test")

ax = nexttile(1);
hold(ax, "on")
colormap('gray')
imagesc(reshape(m_test, 3, 3).')
clim([-0.1 1])
ax.XTick = 1 : 3;
ax.YTick = 1 : 3;
ax.YDir = "reverse";
xlabel("j")
ylabel("i")
ax.XLim = [ax.XTick(1) - 0.5, ax.XTick(end) + 0.5];
ax.YLim = [ax.YTick(1) - 0.5, ax.YTick(end) + 0.5];
axis equal
title("True Model")
colorbar(ax, "eastoutside")

ax = nexttile(2);
hold(ax, "on")
colormap('gray')
imagesc(reshape(m_dagger, 3, 3).')
clim([-0.1 1])
ax.XTick = 1 : 3;
ax.YTick = 1 : 3;
ax.YDir = "reverse";
xlabel("j")
ylabel("i")
ax.XLim = [ax.XTick(1) - 0.5, ax.XTick(end) + 0.5];
ax.YLim = [ax.YTick(1) - 0.5, ax.YTick(end) + 0.5];
axis equal
title("Recoverd Model")
colorbar(ax, "eastoutside")


%% Homework 3, Problem 1
%
% Exercise 2 in Section 3.6

fprintf("**************************\n")
disp("Homework 3, Problem 1")
fprintf("**************************\n\n")

% Set Up Checkerboard Test
m_test = ones(n, 1);
m_test(1:2:n) = -1;
d_test = G * m_test;
disp("Checkboard Test Model")
disp(m_test)
disp("Predicted Data")
disp(d_test)

% Conduct Spike Test
m_dagger = pinv(G) * d_test;
disp("Recovered Model from Checkerboard Test")
disp(m_dagger)

% Vizualize
fig = figure("Name", "Checkerboard Test");
tl = tiledlayout(1, 2, "Parent", fig);
title(tl, "Checkerboard Test")

ax = nexttile(1);
hold(ax, "on")
colormap('gray')
imagesc(reshape(m_test, 3, 3).')
clim([-0.1 1])
ax.XTick = 1 : 3;
ax.YTick = 1 : 3;
ax.YDir = "reverse";
xlabel("j")
ylabel("i")
ax.XLim = [ax.XTick(1) - 0.5, ax.XTick(end) + 0.5];
ax.YLim = [ax.YTick(1) - 0.5, ax.YTick(end) + 0.5];
axis equal
title("True Model")
colorbar(ax, "eastoutside")

ax = nexttile(2);
hold(ax, "on")
colormap('gray')
imagesc(reshape(m_dagger, 3, 3).')
clim([-0.1 1])
ax.XTick = 1 : 3;
ax.YTick = 1 : 3;
ax.YDir = "reverse";
xlabel("j")
ylabel("i")
ax.XLim = [ax.XTick(1) - 0.5, ax.XTick(end) + 0.5];
ax.YLim = [ax.YTick(1) - 0.5, ax.YTick(end) + 0.5];
axis equal
title("Recoverd Model")
colorbar(ax, "eastoutside")

% Examine Model Error
m_error = m_dagger - m_test;
disp("Model Error")
disp(m_error)