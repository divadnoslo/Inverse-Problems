close all
clear
clc

% Save figures as *.eps
saveFigureAsEps = @(name, fig)(exportgraphics(fig, fullfile("..", "latex", "images", name)));

% Seed RNG
rng(0)

%% Task 01: Unknown Polynomial Fit


%% Create Data to Estimate

% Sample Temperatures
temperatureRange = [-50 100];
dT = 0.1;
T = temperatureRange(1) : dT : temperatureRange(end);
ambientTemperature = 21;
sampleTemperatures = -15 : 15 : 60;

% Design Bias Temperature Sensitivity
m_true = [0, 2.1e-4, -3.7e-5, 3.5e-6]; %, 8.4e-8, -9.2e-9];
biasSensitivity = polynomialTemperatureSensitivity(m_true, T, ambientTemperature);

% Create Data
trueSigma = 1e-2;
biasDataTruth = polynomialTemperatureSensitivity(m_true, sampleTemperatures, ambientTemperature);
biasData = biasDataTruth + trueSigma.*randn(size(biasDataTruth));

% Vizualize
fig = figure("Name", "Bias Temperature Sensitivity");
ax = gca;
hold(ax, "on")
plot(T, biasSensitivity, 'k', 'LineWidth', 0.5)
plot(sampleTemperatures, biasData, 'bo')
title("Bias Sensitivity Curve")
xlabel("Temperature [deg C]")
ylabel("[m/sec/sec]")
xlim([temperatureRange(1) temperatureRange(2)])
grid on
grid minor
legend(["True Curve", "Data"], "Location", "Best")


%% Cycle Through Polynomial Fits

% Basic Least Squares Fit
m = length(biasData);
modelOperators = cell(1, m);
models = cell(1, m);
residuals = cell(1, m);
residualNorms = zeros(1, m);
sigmas = zeros(1, m);
conditionNumbers = zeros(1, m);
for n = 1 : m

    % Build Model Operator
    G = zeros(m, n);
    for ii = 1 : m
        for jj = 1 : n
            G(ii,jj) = (sampleTemperatures(ii) - ambientTemperature)^(jj-1);
        end
    end
    modelOperators{n} = G;

    % Perform Least Squares Fit
    m_test = (G.' * G) \ G.' * biasData.';
    models{n} = m_test;

    % Residuals
    residuals{n} = G*m_test - biasData.';
    residualNorms(n) = norm(residuals{n});

    % Sigmas
    sigmas(n) = std(residuals{n});

    % Condition Number
    conditionNumbers(n) = cond(G);

end


%% Analyze Results

colors = ["#8B0000", "#C11C84", "#FF8C00", "#006400", "#00008B", "#301934"];

% Models
fig = figure("Name", "Bias Temperature Sensitivity Models");
ax = gca;
hold(ax, "on")
plot(T, biasSensitivity, 'k', 'LineWidth', 0.5)
plot(sampleTemperatures, biasData, 'bo')
for n = 1 : m
    plot(T, polynomialTemperatureSensitivity(models{n}, T, ambientTemperature), 'Color', colors(n))
end
title("Estimated Models")
xlabel("Temperature [deg C]")
ylabel("[m/sec/sec]")
xlim([sampleTemperatures(1) - 10, sampleTemperatures(end) + 10])
grid on
grid minor
legend(["Truth", "Data", "n = 1", "n = 2", "n = 3", "n = 4", "n = 5", "n = 6",], "Location", "Best")

% Residual Norm
fig = figure("Name", "Residual Norms Comparison");
ax = gca;
hold(ax, "on")
bar(0:m-1, residualNorms)
title("Residual Norms")
xlabel("Polynomial Order")
ylabel("||Gm - d||_2")
xlim([-1, m])
grid on
grid minor

% Sigmas
fig = figure("Name", "Emperical Standard Deviation Comparison");
ax = gca;
hold(ax, "on")
bar(0:m-1, sigmas)
line([0, m+1], [trueSigma trueSigma], 'Color', 'k')
title("Emprical Standard Deviations")
xlabel("Polynomial Order")
ylabel("\sigma")
xlim([-1, m])
grid on
grid minor

% Condition Numbers
fig = figure("Name", "Condition Number Comparison");
ax = gca;
hold(ax, "on")
bar(0:m-1, conditionNumbers)
title("Condition Numbers")
xlabel("Polynomial Order")
ylabel("cond(G)")
xlim([-1, m])
ax.YScale = "log";
grid on
grid minor

