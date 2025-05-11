function out = evalGyro(w_b__i_b_true, w_b__i_b_meas, m_gyro_true, profileNumber, gyroTable, gyroModelLabels)
% Evaluate Gyroscope
%
%

%% Input and Output Arguments

arguments (Input)
    w_b__i_b_true (3,:) {isfloat, mustBeReal, mustBeFinite}
    w_b__i_b_meas (3,:) {isfloat, mustBeReal, mustBeFinite}
    m_gyro_true  (12,1) {isfloat, mustBeReal, mustBeFinite}
    profileNumber (1,1) {mustBePositive}
    gyroTable     (:,:) {mustBeA(gyroTable, 'table')}
    gyroModelLabels
end

arguments (Output)
    out
end


%% Implementation

% Helper Functions
saveFigureAsEps = @(name, fig)(exportgraphics(fig, fullfile("..", "latex", "images", name)));
preamble = @(profileNumber, description)(sprintf("Motion Profile %s: %s", num2str(profileNumber), description));
makeFileName = @(profileNumber, description)(sprintf("MP%s_%s", num2str(profileNumber), description));

% Create Model Operator
Omega = [ones(length(w_b__i_b_true), 1), w_b__i_b_true.'];
Z = zeros(size(Omega));
G = [Omega, Z, Z; Z, Omega, Z; Z, Z, Omega];
[m, n] = size(G);

% Create Data
d = [...
    w_b__i_b_meas(1,:).' - w_b__i_b_true(1,:).'; ...
    w_b__i_b_meas(2,:).' - w_b__i_b_true(2,:).'; ...
    w_b__i_b_meas(3,:).' - w_b__i_b_true(3,:).'];

% Examine Singular Values
[U, S, V] = svd(G, 'econ'); 
fig = figure("Name", "SVD Singular Values");
ax = gca;
hold(ax, "on")
plot(1:n, diag(S), 'bo')
title(preamble(profileNumber, "Gyroscope Singular Values"))
xlabel("Singular Value Index")
ylabel("s_i")
ax.YLim(1) = 1;
ax.YScale = "log";
grid on
grid minor
saveFigureAsEps(makeFileName(profileNumber, "gyro_singular_values.eps"), fig)

% Estimate Calibration Parameters via Normal Equations
m_gyro = inv(G.' * G) * G.' * d;
gyroTable.L2Model = m_gyro;

% Examine Error with True Model Parameters
m_gyro_error = m_gyro - m_gyro_true;
gyroTable.ModelError = m_gyro_error;

% Compute Residual Error
r = d - (G * m_gyro);

% Estimate Standard Deviation
nu = m - n;
s = norm(r) / sqrt(nu);

% Estimate Model Covariance
C_gyro = s^2 * inv(G.' * G);

% Save Variance Values
gyroTable.Covariance = diag(C_gyro);

% Compute 95% Confidence Bound
conf95 = tcdf(0.95, nu) * sqrt(diag(C_gyro));
gyroTable.Confidence95 = conf95;

% Model Correlation Matrix
Rho = computeModelCorrelationMatrix(C_gyro);

% Plot Model Covariance Matrix
fig = figure("Name", "Gyro Correlation matrix");
ax = gca;
colormap('gray')
imagesc(Rho)
clim([-1 1])
ax.XTick = 1 : n;
ax.YTick = 1 : n;
ax.XTickLabel = gyroModelLabels;
ax.YTickLabel = gyroModelLabels;
ax.YDir = "reverse";
axis equal
ax.XLim = [ax.XTick(1) - 0.5, ax.XTick(end) + 0.5];
ax.YLim = [ax.YTick(1) - 0.5, ax.YTick(end) + 0.5];
title(preamble(profileNumber, "Gyro Correlation Matrix"))
colorbar(ax, "eastoutside")
saveFigureAsEps(makeFileName(profileNumber, "gyro_correlation_matrix.eps"), fig)

% Return Results
out = struct();
out.results = gyroTable;
out.singularValues = diag(S);
out.covariance = C_gyro;
out.correlationMatrix = Rho;

% Print Results
fprintf("\n\n")
fprintf("Motion Profile %d: Gyroscope Calibration Results\n", profileNumber)
disp(gyroTable)
fprintf("\n\n")


end

