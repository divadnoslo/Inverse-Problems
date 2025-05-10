%% Task 06: Run Multi-Axis Rate Table Experiment

% close all
% clear
% clc

% Save figures as *.eps
saveFigureAsEps = @(name, fig)(exportgraphics(fig, fullfile("..", "latex", "images", name)));

% Load Working File
load(fullfile(pwd, "working_file.mat"))

% Preamble
preamble = @(description)(sprintf("Motion Profile 3: %s", description));

% Make File Name
makeFileName = @(description)(sprintf("MP3_%s", description));


%% Create Multi-Axis Rate Table Motion

g = [0; 0; -gravity];

Fs = 100;
dt = 1 / Fs;
t = 0 : dt : 10;
K = length(t);

w_b__i_b_true = zeros(3, K);
w_b__i_b_true(1,:) = 30 * (2*pi) * (pi/180) * sin(2*pi*t);
w_b__i_b_true(2,:) = 45 * (2*pi) * (pi/180) * cos(2*2*pi*t);
w_b__i_b_true(3,:) = 60 * (2*pi) * (pi/180) * sin(pi*t);

euler = cumtrapz(t, w_b__i_b_true, 2);
dcm = euler2dcm(euler(1,:), euler(2,:), euler(3,:));

f_b__i_b_true = squeeze(pagemtimes(dcm, g));


%% Create IMU Measurements

[f_b__i_b_meas, w_b__i_b_meas] = imu.runForwardModel(f_b__i_b_true, w_b__i_b_true, 21);


%% Vizualize

% Vizualize Angular Velocity
fig = figure("Name", "Angular Velocity Profile");
tl = tiledlayout(3, 1, "Parent", fig);
title(tl, preamble("Angular Velocity Profile"))
ax = nexttile(1);
hold(ax, "on")
plot(t, 180/pi * w_b__i_b_true(1,:), 'k', 'LineWidth', 2)
plot(t, 180/pi * w_b__i_b_meas(1,:), 'r')
title("\omega_x")
xlabel("Time [sec]")
ylabel("[deg/sec]")
xlim([t(1) t(end)])
grid on
grid minor
legend(["Test Bed", "UUT"], "Location", "eastoutside")
ax = nexttile(2);
hold(ax, "on")
plot(t, 180/pi * w_b__i_b_true(2,:), 'k', 'LineWidth', 2)
plot(t, 180/pi * w_b__i_b_meas(2,:), 'g')
title("\omega_y")
xlabel("Time [sec]")
ylabel("[deg/sec]")
xlim([t(1) t(end)])
grid on
grid minor
legend(["Test Bed", "UUT"], "Location", "eastoutside")
ax = nexttile(3);
hold(ax, "on")
plot(t, 180/pi * w_b__i_b_true(3,:), 'k', 'LineWidth', 2)
plot(t, 180/pi * w_b__i_b_meas(3,:), 'b')
title("\omega_z")
xlabel("Time [sec]")
ylabel("[deg/sec]")
xlim([t(1) t(end)])
grid on
grid minor
legend(["Test Bed", "UUT"], "Location", "eastoutside")
linkaxes(tl.Children, 'x')
saveFigureAsEps(makeFileName("angular_velocity_profile.eps"), fig)

% Vizualize Angular Velocity Error
fig = figure("Name", "Angular Velocity Error");
tl = tiledlayout(3, 1, "Parent", fig);
title(tl, preamble("Angular Velocity Error"))
ax = nexttile(1);
hold(ax, "on")
plot(t, 180/pi * (w_b__i_b_meas(1,:) - w_b__i_b_true(1,:)), 'r', 'LineWidth', 2)
title("\Delta\omega_x")
xlabel("Time [sec]")
ylabel("[deg/sec]")
xlim([t(1) t(end)])
grid on
grid minor
ax = nexttile(2);
hold(ax, "on")
plot(t, 180/pi * (w_b__i_b_meas(2,:) - w_b__i_b_true(2,:)), 'r', 'LineWidth', 2)
title("\Delta\omega_y")
xlabel("Time [sec]")
ylabel("[deg/sec]")
xlim([t(1) t(end)])
grid on
grid minor
ax = nexttile(3);
hold(ax, "on")
plot(t, 180/pi * (w_b__i_b_meas(3,:) - w_b__i_b_true(3,:)), 'r', 'LineWidth', 2)
title("\Delta\omega_z")
xlabel("Time [sec]")
ylabel("[deg/sec]")
xlim([t(1) t(end)])
grid on
grid minor
linkaxes(tl.Children, 'x')
saveFigureAsEps(makeFileName("angular_velocity_error.eps"), fig)

% Vizualize Euler Angle Profile
fig = figure("Name", "Euler Angle Profile");
tl = tiledlayout(3, 1, "Parent", fig);
title(tl, preamble("Euler Angle Profile"))
ax = nexttile(1);
hold(ax, "on")
plot(t, 180/pi * euler(1,:), 'r')
title("Roll (\phi)")
xlabel("Time [sec]")
ylabel("[deg]")
xlim([t(1) t(end)])
grid on
grid minor
ax = nexttile(2);
hold(ax, "on")
plot(t, 180/pi * euler(2,:), 'g')
title("Pitch (\theta)")
xlabel("Time [sec]")
ylabel("[deg]")
xlim([t(1) t(end)])
grid on
grid minor
ax = nexttile(3);
hold(ax, "on")
plot(t, 180/pi * euler(3,:), 'b')
title("Yaw (\psi)")
xlabel("Time [sec]")
ylabel("[deg]")
xlim([t(1) t(end)])
grid on
grid minor
linkaxes(tl.Children, 'x')
saveFigureAsEps(makeFileName("euler_angle_profile.eps"), fig)

% Vizualize Specific Force Profile
fig = figure("Name", "Specific Force Profile");
tl = tiledlayout(3, 1, "Parent", fig);
title(tl, preamble("Specific Force Profile"))
ax = nexttile(1);
hold(ax, "on")
plot(t, f_b__i_b_true(1,:), 'k', 'LineWidth', 2)
plot(t, f_b__i_b_meas(1,:), 'r')
title("f_x")
xlabel("Time [sec]")
ylabel("[m/sec/sec]")
xlim([t(1) t(end)])
grid on
grid minor
legend(["Test Bed", "UUT"], "Location", "eastoutside")
ax = nexttile(2);
hold(ax, "on")
plot(t, f_b__i_b_true(2,:), 'k', 'LineWidth', 2)
plot(t, f_b__i_b_meas(2,:), 'g')
title("f_y")
xlabel("Time [sec]")
ylabel("[m/sec/sec]")
xlim([t(1) t(end)])
grid on
grid minor
legend(["Test Bed", "UUT"], "Location", "eastoutside")
ax = nexttile(3);
hold(ax, "on")
plot(t, f_b__i_b_true(3,:), 'k', 'LineWidth', 2)
plot(t, f_b__i_b_meas(3,:), 'b')
title("f_z")
xlabel("Time [sec]")
ylabel("[m/sec/sec]")
xlim([t(1) t(end)])
grid on
grid minor
legend(["Test Bed", "UUT"], "Location", "eastoutside")
linkaxes(tl.Children, 'x')
saveFigureAsEps(makeFileName("specific_force_profile.eps"), fig)

% Vizualize Specific Force Error
fig = figure("Name", "Specific Force Error");
tl = tiledlayout(3, 1, "Parent", fig);
title(tl, preamble("Specific Force Error"))
ax = nexttile(1);
hold(ax, "on")
plot(t, f_b__i_b_meas(1,:) - f_b__i_b_true(1,:), 'r', 'LineWidth', 2)
title("\Delta f_x")
xlabel("Time [sec]")
ylabel("[m/sec/sec]")
xlim([t(1) t(end)])
grid on
grid minor
ax = nexttile(2);
hold(ax, "on")
plot(t, f_b__i_b_meas(2,:) - f_b__i_b_true(2,:), 'r', 'LineWidth', 2)
title("\Delta f_y")
xlabel("Time [sec]")
ylabel("[m/sec/sec]")
xlim([t(1) t(end)])
grid on
grid minor
ax = nexttile(3);
hold(ax, "on")
plot(t, f_b__i_b_meas(3,:) - f_b__i_b_true(3,:), 'r', 'LineWidth', 2)
title("\Delta f_z")
xlabel("Time [sec]")
ylabel("[m/sec/sec]")
xlim([t(1) t(end)])
grid on
grid minor
linkaxes(tl.Children, 'x')
saveFigureAsEps(makeFileName("specific_force_error.eps"), fig)


%% Evaluate

mp3Gyro = evalGyro(...
    w_b__i_b_true, ...
    w_b__i_b_meas, ...
    m_gyro_true, ...
    3, ...
    gyroTable, ...
    gyroModelLabels)

mp3Accel = evalAccel(...
    f_b__i_b_true, ...
    f_b__i_b_meas, ...
    m_accel_true, ...
    3, ...
    accelTable, ...
    accelModelLabels)

%% Append to Working File

save(workingFilePath, "w_b__i_b_true", "f_b__i_b_true", "w_b__i_b_meas", "f_b__i_b_meas", "-append")

