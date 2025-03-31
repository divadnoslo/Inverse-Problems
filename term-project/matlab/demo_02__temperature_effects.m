close all
clear
clc

% Save figures as *.eps
saveFigureAsEps = @(name, fig)(exportgraphics(fig, fullfile("..", "latex", "images", name)));


%% Demo 02: Temperature Effects

% Sample Temperatures
temperatureRange = [-50 100];
dT = 0.1;
T = temperatureRange(1) : dT : temperatureRange(end);
ambientTemperature = 21;

% Design Bias Temperature Sensitivity
m_true = [0, 3.4e-5, 7.5e-6, 2.4e-7];
biasSensitivity = polynomialTemperatureSensitivity(m_true, T, ambientTemperature);

fig = figure("Name", "Bias Temperature Sensitivity");
plot(T, biasSensitivity, 'b')
title("Bias Sensitivity Curve")
xlabel("Temperature [deg C]")
ylabel("[m/sec/sec]")
grid on
grid minor

% Simulate an IMU with Bias Temperature Sensitivity
imu = ImuModel;

% Set Accelerometer Parameters
imu.b_a = [0.01; 0.02; 0.03];
imu.M_a = [...
    0.001,   0.0001,  0.0002; ...
    0.0003,  0.002,   0.0004; ...
    0.0005,  0.0006,  0.003]; 
imu.sigma_a = [3e-2, 4e-2, 5e-2];

% Set Gyroscope Parameters
imu.b_g = [0.01; 0.02; 0.03];
imu.M_g = [...
    0.001,   0.0001,  0.0002; ...
    0.0003,  0.002,   0.0004; ...
    0.0005,  0.0006,  0.003];  
imu.sigma_g = [3e-2, 4e-2, 5e-2] .* pi/180;

% Set Accel Bias Sensitivity
imu.AccelBiasTempSensitivityModel = @(temp)(polynomialTemperatureSensitivity(m_true, temp, ambientTemperature));


%% Simulate Multiple Calibration Cycles

icm = ImuCalibrationManager(...
    "SamplingRate", 100, ...
    "Duration", 10);

sampleTemps = temperatureRange(1) : 5 : temperatureRange(end);
biases = zeros(1, length(sampleTemps));

for k = 1 : length(biases)

    dataset = icm.createCalibrationDataset(imu, sampleTemps(k));
    newImu = icm.processCalibrationDataset(dataset);
    biases(k) = newImu.AccelFixedBiasX;

end

fig = figure("Name", "Example Bias Sensitivity Characterization");
plot(sampleTemps, biases, 'bo')
title("Example Bias Sensitivity Characterization")
xlabel("Temperature [deg C]")
ylabel("[m/sec/sec]")
grid on
grid minor
saveFigureAsEps("ideal_bias_response_to_temperature.eps", fig)

