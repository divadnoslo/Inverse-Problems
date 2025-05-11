%% Task 02: Run Single-Axis Rate Table Experiment

% close all
% clear
% clc

% Load Working File
load(fullfile(pwd, "working_file.mat"))

% Profile Number
profileNumber = 2;

% Helper Functions
saveFigureAsEps = @(name, fig)(exportgraphics(fig, fullfile("..", "latex", "images", name)));
preamble = @(profileNumber, description)(sprintf("Motion Profile %s: %s", num2str(profileNumber), description));
makeFileName = @(profileNumber, description)(sprintf("MP%s_%s", num2str(profileNumber), description));


%% Create Single-Axis Rate-Table Motion

g = [0; 0; -gravity];

Fs = 100;
dt = 1 / Fs;
t_end = 5;
t_segment = 0 : dt : t_end;
K = length(t_segment);

A = (2*pi) * 45 * pi/180;
curve = A * sin(2*pi*(1/t_end)*t_segment);

curve = [curve, -curve];
K = length(curve);

w_b__i_b_true = zeros(3, 3*K);
w_b__i_b_true(1,1:K) = curve;
w_b__i_b_true(2,K + (1:K)) = curve;
w_b__i_b_true(3,2*K + (1:K)) = curve;

t = (0:length(w_b__i_b_true)-1) * dt;

euler = cumtrapz(t, w_b__i_b_true, 2);
dcm = euler2dcm(euler(1,:), euler(2,:), euler(3,:));

f_b__i_b_true = squeeze(pagemtimes(dcm, g));


%% Create IMU Measurements

[f_b__i_b_meas, w_b__i_b_meas] = imu.runForwardModel(f_b__i_b_true, w_b__i_b_true, 21);


%% Vizualize

% Vizualize Angular Velocity
fig = figure("Name", "Angular Velocity Profile");
tl = tiledlayout(3, 1, "Parent", fig);
title(tl, preamble(profileNumber, "Angular Velocity Profile"))
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
saveFigureAsEps(makeFileName(profileNumber, "angular_velocity_profile.eps"), fig)

% Vizualize Angular Velocity Error
fig = figure("Name", "Angular Velocity Error");
tl = tiledlayout(3, 1, "Parent", fig);
title(tl, preamble(profileNumber, "Angular Velocity Error"))
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
saveFigureAsEps(makeFileName(profileNumber, "angular_velocity_error.eps"), fig)

% Vizualize Euler Angle Profile
fig = figure("Name", "Euler Angle Profile");
tl = tiledlayout(3, 1, "Parent", fig);
title(tl, preamble(profileNumber, "Euler Angle Profile"))
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
saveFigureAsEps(makeFileName(profileNumber, "euler_angle_profile.eps"), fig)

% Vizualize Specific Force Profile
fig = figure("Name", "Specific Force Profile");
tl = tiledlayout(3, 1, "Parent", fig);
title(tl, preamble(profileNumber, "Specific Force Profile"))
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
saveFigureAsEps(makeFileName(profileNumber, "specific_force_profile.eps"), fig)

% Vizualize Specific Force Error
fig = figure("Name", "Specific Force Error");
tl = tiledlayout(3, 1, "Parent", fig);
title(tl, preamble(profileNumber, "Specific Force Error"))
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
saveFigureAsEps(makeFileName(profileNumber, "specific_force_error.eps"), fig)


%% Evaluate

mp2Gyro = evalGyro(...
    w_b__i_b_true, ...
    w_b__i_b_meas, ...
    m_gyro_true, ...
    profileNumber, ...
    gyroTable, ...
    gyroModelLabels);

mp2Accel = evalAccel(...
    f_b__i_b_true, ...
    f_b__i_b_meas, ...
    m_accel_true, ...
    profileNumber, ...
    accelTable, ...
    accelModelLabels);

%% Append to Working File

save(workingFilePath, "mp2Gyro", "mp2Accel", "-append")

