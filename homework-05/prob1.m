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

fh = @(m, t)(m(1)*sin(2*pi*m(2).*t + m(3)) + m(4));

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
A0 = (max(instdata - c0) - min(instdata - c0)) / 4;
chi20 = -11*pi/8;

m0 = [A0; f0; chi20; c0];

% Print Initial Model Guess
fprintf("Initial Model Guess\n")
fprintf("\tA0:    %6.3f V\n",  A0)
fprintf("\tf0:    %6.3f Hz\n", f0)
fprintf("\tchi20: %6.3f   \n",  chi20)
fprintf("\tc0:    %6.3f V\n",  c0)

% Vizually Inspect Initial Guess
fig = figure("Name", "Initial Model Guess");
ax = gca;
hold(ax, "on")
plot(t, instdata, 'b.')
plot(t, fh(m0,t), 'r', 'LineWidth', 3)
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
fprintf("\tchi2:  %6.3f  \n",  m(3))
fprintf("\tc:     %6.3f V\n",  m(4))

% Vizually Inspect Levenberg-Marquardt Solution
fig = figure("Name", "Levenberg-Marquardt Solution");
ax = gca;
hold(ax, "on")
plot(t, instdata, 'b.')
plot(t, fh(m,t), 'm', 'LineWidth', 3)
title("Levenberg-Marquardt Solution")
xlabel("Time [sec]")
ylabel("[V]")
grid on
grid minor
saveFigureAsEps("prob1_lm_solution.eps", fig)

% Estimate Additive Noise
residuals = instdata - prob1Function(m);
s = sqrt((norm(residuals, 2).^2 ./ (length(instdata) - length(m))));
fprintf("\nEstimated Standard Deviation\n")
fprintf("\ts = %6.3f\n", s)

% Plot Residuals
fig = figure("Name", "Levenberg-Marquardt Residuals");
ax = gca;
hold(ax, "on")
p = patch([0 40 40 0], [-3*s, -3*s, 3*s, 3*s], 'c');
p.FaceAlpha = 0.2;
p.EdgeColor = 'none';
plot(t, residuals, 'b.')
line([ax.XLim(1) ax.XLim(2)], [mean(residuals) mean(residuals)], 'Color', 'k', 'LineStyle', '--')
title("Levenberg-Marquardt Residuals")
xlabel("Time [sec]")
ylabel("[V]")
grid on
grid minor
text(5, 3*s - 0.1*3*s, sprintf("3-Sigma Bound = %6.3f", 3*s), "FontWeight", "bold")
saveFigureAsEps("prob1_lm_solution_residuals.eps", fig)

% Plot Historgam of Residuals for Kicks
fig = figure("Name", "Histogram of LM Residuals");
histogram(residuals, "Normalization", "pdf", "BinWidth", 25)
title("Histogram of LM Residuals")
xlabel("Residuals")
ylabel("pdf")
grid on
grid minor
saveFigureAsEps("prob1_lm_residual_histogram.eps", fig)

% Model Covariance Estimate
J = prob1Jacobian(m);
C = s^2 * inv(J.' * J);
fprintf("\nEstimated Model Covariance\n")
disp(C)

% 95% Confidence Interval
cf95 = 1.96 .* sqrt(diag(C));
fprintf("95%% Confidence Intervals\n")
fprintf("\tA:    %6.3f +/- %6.6f V\n", m(1), cf95(1))
fprintf("\tf0:   %6.3f +/- %6.6f V\n", m(2), cf95(2))
fprintf("\tchi2: %6.3f +/- %6.6f  \n", m(3), cf95(3))
fprintf("\tc:    %6.3f +/- %6.6f V\n", m(4), cf95(4))

% Model Paramter Correlation Matrix
Rho = zeros(4, 4);
for ii = 1 : 4
    for jj = 1 : 4
        Rho(ii,jj) = C(ii,jj) ./ sqrt(C(ii,ii) * C(jj,jj));
    end
end
fprintf("\nCorrelation Matrix\n")
disp(Rho)


%% Hypothetical Bad Starting Model

% Perform Levenberg-Marqaurdt Estimation
[m_bad, iter] = lm('prob1Function', 'prob1Jacobian', zeros(4, 1), 1e-18, 100);

% Print Levenberg-Marquardt Solution
fprintf("\nBad Levenberg-Marquardt Solution\n")
fprintf("\tA:     %6.3f V\n",  m_bad(1))
fprintf("\tf0:    %6.3f Hz\n", m_bad(2))
fprintf("\tchi2:  %6.3f  \n",  m_bad(3))
fprintf("\tc:     %6.3f V\n",  m_bad(4))

% Vizually Inspect Levenberg-Marquardt Solution
fig = figure("Name", "Levenberg-Marquardt Solution");
ax = gca;
hold(ax, "on")
plot(t, instdata, 'b.')
plot(t, fh(m_bad, t), 'm', 'LineWidth', 3)
title("Levenberg-Marquardt Solution (Bad Starting Model)")
xlabel("Time [sec]")
ylabel("[V]")
grid on
grid minor
saveFigureAsEps("prob1_lm_solution_bad_start_model.eps", fig)


%% Other Good Models

m_alt = [10; 0.1; -25; -18];

% Perform Levenberg-Marqaurdt Estimation
[m_new, iter] = lm('prob1Function', 'prob1Jacobian', m_alt, 1e-24, 100);

% Print Levenberg-Marquardt Solution
fprintf("\nAlternative Levenberg-Marquardt Solution\n")
fprintf("\tA:     %6.3f V\n",  m_new(1))
fprintf("\tf0:    %6.3f Hz\n", m_new(2))
fprintf("\tchi2:  %6.3f  \n",  m_new(3))
fprintf("\tc:     %6.3f V\n",  m_new(4))

% Vizually Inspect Levenberg-Marquardt Solution
fig = figure("Name", "Levenberg-Marquardt Solution");
ax = gca;
hold(ax, "on")
plot(t, instdata, 'b.')
plot(t, fh(m_alt, t), 'r', 'LineWidth', 3)
plot(t, fh(m_new, t), 'm', 'LineWidth', 3)
title("Levenberg-Marquardt Solution (Alternative Starting Model)")
xlabel("Time [sec]")
ylabel("[V]")
grid on
grid minor
legend(["Initial Model", "LM Solution"], "Location", "Best")
saveFigureAsEps("prob1_lm_solution_new_start_model.eps", fig)


