close all
clear
clc

%% Lecture Example on 24 Jan 2025

% A ball is thrown upward with an initial velocity and comes back down
m_true = [...
    10; ...  % initial height [m]
    3; ...   % initial velocity [m/sec]
    9.81];   % acceleration due to gravity [m/sec/sec]

trueModel = @(t, m)(m(1) + m(2).*t - (1/2).*m(3).*t.^2);

% Create observations with noise
t = 1 : 10;
t = t.';
d_true = trueModel(t, m_true);
oneSigma = 5e-2;
d_meas = d_true + oneSigma .* randn(size(d_true));

% Create G
G = zeros(10, 3);
for k = 1 : 10
    G(k,:) = [1, t(k), -(1/2)*t(k)^2];
end

% Solve via Least Squares
m_est = (G.' * G) \ (G.' * d_meas);

% Print Results
fprintf("Initial Height\n")
fprintf("\tTruth:    %8.6f m\n", m_true(1));
fprintf("\tEstimate: %8.6f m\n", m_est(1));
fprintf("\tError:    %8.3e m\n", m_est(1) - m_true(1))
fprintf("\n")
fprintf("Initial Velocity\n")
fprintf("\tTruth:    %8.6f m/sec\n", m_true(2));
fprintf("\tEstimate: %8.6f m/sec\n", m_est(2));
fprintf("\tError:    %8.3e m/sec\n", m_est(2) - m_true(2))
fprintf("\n")
fprintf("Acceleration due to Gravity\n")
fprintf("\tTruth:    %8.6f m/sec/sec\n", m_true(3));
fprintf("\tEstimate: %8.6f m/sec/sec\n", m_est(3));
fprintf("\tError:    %8.3e m/sec/sec\n", m_est(3) - m_true(3))
fprintf("\n")

