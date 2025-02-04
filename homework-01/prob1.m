close all
clear
clc

% Problem 1

%% Part A

% Discritize the integral into 20 equally spaced points
n = 20;
dz = 1/n;
z = zeros(n, 1);
for k = 1 : n
    z(k) = dz/2 + (k - 1)*dz;
end
x = z;

% Function handles to "g_{i,j}" and "d_{i}"
g_ij = @(x, z)(5*sin(x.*z));
d_i = @(x)(50*sin(x) - 50*sin(x).*cos(x));

% Construct "G" and "d" quantities
G = zeros(n, n);
for i = 1 : n
    for j = 1 : n
        G(i,j) = g_ij(x(i), z(j));
    end
end

d = zeros(n, 1);
for i = 1 : 20
    d(i) = d_i(x(i));
end

% Solve for "m" in "Gm = d"
m = G \ d;

% Print results to the command window
fprintf("**********************\n")
fprintf("Part A\n")
fprintf("**********************\n\n")
fprintf("G Matrix\n")
disp(G)
fprintf("\n\nData (d)\n")
disp(d)
fprintf("\n\nModel (m)\n")
disp(m)


%% Part B

% Compute condition number of G
conditionNumber = cond(G);

% Print condition number
fprintf("**********************\n")
fprintf("Part B\n")
fprintf("**********************\n\n")
fprintf("Condition number of G: %10.6e\n\n", conditionNumber)


%% Part C

% True Model
trueModel = @(x)(10.*x.*sin(x));
m_true = zeros(n, 1);
for i = 1 : n
    m_true(i) = trueModel(x(i));
end

% Plot Results
fig = figure("Name", "True vs. Computed Model");
tl = tiledlayout(2, 1, "Parent", fig);

nexttile(1)
hold on
plot(x, m_true, 'k.', "MarkerSize", 5)
plot(x, m, 'b.', 'MarkerSize', 5)
title("True vs. Computed Model")
xlabel("x")
ylabel("m")
grid on
grid minor
legend(["True Model", "Computed Model"], "Location", "Best")
hold off

nexttile(2)
hold on
plot(x, m - m_true, 'r.')
title("Error in Computed Model")
xlabel("x")
ylabel("\Delta m")
grid on
grid minor
hold off

linkaxes(tl.Children, 'x')