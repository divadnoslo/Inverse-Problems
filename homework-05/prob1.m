close all
clear
clc

% Add to Path
addpath(genpath(fullfile("..", "PEIP-master", "Lib")))

% Save figures as *.eps
saveFigureAsEps = @(name, fig)(exportgraphics(fig, fullfile(pwd, "latex", "images", name)));

% Load Data
load(fullfile(pwd, "data", "instdata.mat"))


%% Problem 1 - Exercise 2 in Section 9.6

n = 1 : length(instdata);
Fs = 50;
dt = 1 / Fs;
t = (n - 1) * dt;

% Plot Data
fig = figure("Name", "Instrument Recording");
plot(t, instdata, 'b.')
title("Provided Instrument Recording")
xlabel("Time [sec]")
ylabel("[V]")
grid on
grid minor
saveFigureAsEps("prob1_instrument_recording.eps", fig)

% Plot PSD for Kicks
[pxx, f] = pwelch(instdata, [], [], 2^10, Fs, "onesided");
fig = figure("Name", "PSD of Instrument Data");
plot(f, 10*log10(pxx))
title("PSD of Instrument Recording")
xlabel("Frequency [Hz]")
ylabel("[dB]")
grid on
grid minor
saveFigureAsEps("prob1_instrument_recording_PSD.eps", fig)

% Take a Guess at an Initial Model
[~, psdMaxIndex] = max(pxx);
f0 = f(psdMaxIndex);
c0 = mean(instdata);
A0 = (max(instdata - c0) - min(instdata - c0)) / 2;
chi20 = pi/36;

m0 = [A0; f0; chi20; c0];

% Print Initial Model Guess
fprintf("Initial Model Guess\n")
fprintf("\tA0:    %6.3f V\n",  A0)
fprintf("\tf0:    %6.3f Hz\n", f0)
fprintf("\tchi20: %6.3f V\n",  chi20)
fprintf("\tc0:    %6.3f V\n",  c0)

% Vizually Inspect Initial Guess
fig = figure("Name", "Initial Model Guess");
ax = gca;
hold(ax, "on")
plot(t, instdata, 'b.')
plot(t, prob1Function(m0), 'r', 'LineWidth', 3)
title("Vizually Inspecting Initial Model Guess")
xlabel("Time [sec]")
ylabel("[V]")
grid on
grid minor
saveFigureAsEps("prob1_initial_model_guess.eps", fig)

% Perform Levenberg-Marqaurdt Estimation
[m, iter] = lm('prob1Function', 'prob1Jacobian', m0, 1e-18, 100);

% Print Levenberg-Marquardt Solution
fprintf("\nLevenberg-Marquardt Solution\n")
fprintf("\tA:     %6.3f V\n",  m(1))
fprintf("\tf0:    %6.3f Hz\n", m(2))
fprintf("\tchi2:  %6.3f V\n",  m(3))
fprintf("\tc:     %6.3f V\n",  m(4))

% Vizually Inspect Levenberg-Marquardt Solution
fig = figure("Name", "Levenberg-Marquardt Solution");
ax = gca;
hold(ax, "on")
plot(t, instdata, 'b.')
plot(t, prob1Function(m), 'm', 'LineWidth', 3)
title("Levenberg-Marquardt Solution")
xlabel("Time [sec]")
ylabel("[V]")
grid on
grid minor
saveFigureAsEps("prob1_lm_solution.eps", fig)

% Do This While LM is Broken
m = m0;

% Estimate Additive Noise
residuals = instdata - prob1Function(m);
s = std(residuals);

% Plot Residuals
fig = figure("Name", "Levenberg-Marquardt Residuals");
ax = gca;
hold(ax, "on")
p = patch([0 40 40 0], [-3*s, -3*s, 3*s, 3*s], 'c');
p.FaceAlpha = 0.2;
p.EdgeColor = 'none';
plot(t, residuals, 'b.')
title("Levenberg-Marquardt Residuals")
xlabel("Time [sec]")
ylabel("[V]")
grid on
grid minor
text(5, 3*s - 0.1*3*s, sprintf("3-Sigma Bound = %6.3f", 3*s), "FontWeight", "bold")
saveFigureAsEps("prob1_lm_solution_residuals.eps", fig)

% Plot PSD of Residuals for Kicks
[rxx, f] = pwelch(residuals, [], [], 2^10, Fs, "onesided");
fig = figure("Name", "PSD of LM Residuals");
plot(f, 10*log10(rxx))
title("PSD of LM Residuals")
xlabel("Frequency [Hz]")
ylabel("[dB]")
grid on
grid minor
saveFigureAsEps("prob1_lm_residual_PSD.eps", fig)

% Model Covariance Estimate
J = prob1Jacobian(m);
C = inv(J.' * J);
fprintf("\nEstimated Model Covariance\n")
disp(C)

% 95% Confidence Interval
cf95 = 1.96 .* sqrt(diag(C));
fprintf("95%% Confidence Intervals\n")
fprintf("\tA:    %6.3f +/- %6.3f V\n", m(1), cf95(1))
fprintf("\tf0:   %6.3f +/- %6.3f V\n", m(2), cf95(2))
fprintf("\tchi2: %6.3f +/- %6.3f V\n", m(3), cf95(3))
fprintf("\tc:    %6.3f +/- %6.3f V\n", m(4), cf95(4))

% Model Paramter Correlation Matrix
Rho = zeros(4, 4);
for ii = 1 : 4
    for jj = 1 : 4
        Rho(ii,jj) = C(ii,jj) ./ sqrt(C(ii,ii) * C(jj,jj));
    end
end
fprintf("\nCorrelation Matrix\n")
disp(Rho)

