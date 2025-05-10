function out = evalAccel(f_b__i_b_true, f_b__i_b_meas, m_accel_true, profileNumber, accelTable, accelModelLabels)
% Evaluate Accelerometers
%
%

%% Input and Output Arguments

arguments (Input)
    f_b__i_b_true (3,:) {isfloat, mustBeReal, mustBeFinite}
    f_b__i_b_meas (3,:) {isfloat, mustBeReal, mustBeFinite}
    m_accel_true  (12,1) {isfloat, mustBeReal, mustBeFinite}
    profileNumber (1,1) {mustBePositive}
    accelTable     (:,:) {mustBeA(accelTable, 'table')}
    accelModelLabels
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
F = [ones(length(f_b__i_b_true), 1), f_b__i_b_true.'];
Z = zeros(size(F));
G = [F, Z, Z; Z, F, Z; Z, Z, F];
[m, n] = size(G);

% Create Data
d = [...
    f_b__i_b_meas(1,:).' - f_b__i_b_true(1,:).'; ...
    f_b__i_b_meas(2,:).' - f_b__i_b_true(2,:).'; ...
    f_b__i_b_meas(3,:).' - f_b__i_b_true(3,:).'];

% Examine Singular Values
[U, S, V] = svd(G, 'econ'); 
fig = figure("Name", "SVD Singular Values");
ax = gca;
hold(ax, "on")
plot(1:n, diag(S), 'bo')
title(preamble(profileNumber, "Accelerometer Singular Values"))
xlabel("Singular Value Index")
ylabel("s_i")
ax.YLim(1) = 0;
ax.YScale = "log";
grid on
grid minor
saveFigureAsEps(makeFileName(profileNumber, "accel_singular_values.eps"), fig)

% Estimate Calibration Parameters via Normal Equations
m_accel = inv(G.' * G) * G.' * d;
accelTable.L2Model = m_accel;

% Examine Error with True Model Parameters
m_accel_error = m_accel - m_accel_true;
accelTable.ModelError = m_accel_error;

% Model Covariance
C_accel = inv(G.' * G);

% Save Variance Values
accelTable.Covariance = diag(C_accel);

% Compute 95% Confidence Bound
conf95 = 1.96 * sqrt(diag(C_accel));
accelTable.Confidence95 = conf95;

% Model Correlation Matrix
Rho = computeModelCorrelationMatrix(C_accel);

% Plot Model Covariance Matrix
fig = figure("Name", "Accel Correlation matrix");
ax = gca;
colormap('gray')
imagesc(Rho)
clim([-1 1])
ax.XTick = 1 : n;
ax.YTick = 1 : n;
ax.XTickLabel = accelModelLabels;
ax.YTickLabel = accelModelLabels;
ax.YDir = "reverse";
axis equal
ax.XLim = [ax.XTick(1) - 0.5, ax.XTick(end) + 0.5];
ax.YLim = [ax.YTick(1) - 0.5, ax.YTick(end) + 0.5];
title(preamble(profileNumber, "Accel Correlation Matrix"))
colorbar(ax, "eastoutside")
saveFigureAsEps(makeFileName(profileNumber, "accel_correlation_matrix.eps"), fig)

% Return Results
out = struct();
out.results = accelTable;
out.singularValues = diag(S);
out.covariance = C_accel;
out.correlationMatrix = Rho;


end

