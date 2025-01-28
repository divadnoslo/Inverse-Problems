close all
clear
clc

%% Problem 01

n = 20;
dx = 1/n;
x = 0 : dx : 1;

kernal = @(x, z)(5.*sin( x .* z));
dataModel = @(x)(50.*sin(x) - 50.*sin(x).*cos(x));
trueModel = @(x)(10.*x.*sin(x));





% Plot Results
fig = figure("Name", "IFK vs. True Model");
hold on
plot(x, trueModel(x), 'k', "LineWidth", 2)
plot(x, dataModel(x), 'r.', 'MarkerSize', 5)
title("IFK vs. True Model")
xlabel("x")
ylabel("m_t_r_u_e(x)")
grid on
grid minor
hold off