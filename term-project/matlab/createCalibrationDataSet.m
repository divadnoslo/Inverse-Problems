function dataset = createCalibrationDataSet(imu, options)
% Create Calibration Data Set
%
%

%% Input and Output Arguments

arguments (Input)
    imu (1,1) ImuInterface
    options.SamplingRate (1,1) {isfloat, mustBeReal, mustBePositive} = 100
    options.Duration (1,1) {isfloat, mustBeReal, mustBePositive} = 10
    options.Gravity (1,1) {isfloat, mustBeReal, mustBePositive} = 9.81
    options.SpinRate (1,1) {isfloat, mustBeReal, mustBeFinite} = 30 * pi/180
    options.Temperature (1,1) {isfloat, mustBeReal, mustBeFinite} = 21
end

arguments (Output)
    dataset (1,1) {mustBeA(dataset, 'struct')}
end


%% Implementation

dataset = struct();

% Test Names
testNames = [...
    "accelPosX", ...
    "accelNegX", ...
    "accelPosY", ...
    "accelNegY", ...
    "accelPosZ", ...
    "accelNegZ", ...
    "gyroPosX", ...
    "gyroNegX", ...
    "gyroPosY", ...
    "gyroNegY", ...
    "gyroPosZ", ...
    "gyroNegZ"];

% Simulate Tests
dataset.time = 0 : (1/options.SamplingRate) : options.Duration;
K = length(dataset.time);
for k = 1 : length(testNames)

    % Initialize
    f_b__i_b_true = zeros(3, K);
    w_b__i_b_true = zeros(3, K);

    % Simulate Perfect Test
    switch testNames(k)

        case "accelPosX"
            f_b__i_b_true(1,:) = options.Gravity * ones(1, K);

        case "accelNegX"
            f_b__i_b_true(1,:) = -1 * options.Gravity * ones(1, K);

        case "accelPosY"
            f_b__i_b_true(2,:) = options.Gravity * ones(1, K);

        case "accelNegY"
            f_b__i_b_true(2,:) = -1 * options.Gravity * ones(1, K);

        case "accelPosZ"
            f_b__i_b_true(3,:) = options.Gravity * ones(1, K);

        case "accelNegZ"
            f_b__i_b_true(3,:) = -1 * options.Gravity * ones(1, K);

        case "gyroPosX"
            w_b__i_b_true(1,:) = options.SpinRate * ones(1, K);

        case "gyroNegX"
            w_b__i_b_true(1,:) = -1 * options.SpinRate * ones(1, K);

        case "gyroPosY"
            w_b__i_b_true(2,:) = options.SpinRate * ones(1, K);

        case "gyroNegY"
            w_b__i_b_true(2,:) = -1 * options.SpinRate * ones(1, K);

        case "gyroPosZ"
            w_b__i_b_true(3,:) = options.SpinRate * ones(1, K);

        case "gyroNegZ"
            w_b__i_b_true(3,:) = -1 * options.SpinRate * ones(1, K);

        otherwise
            error("Test not recognized!")

    end

    % Run IMU Forward Model
    [f_b__i_b_meas, w_b__i_b_meas] = imu.runForwardModel(...
        f_b__i_b_true, ...
        w_b__i_b_true);

    % Store Test Data
    dataset.(testNames(k)) = struct();
    dataset.(testNames(k)).f_b__i_b_true = f_b__i_b_true;
    dataset.(testNames(k)).w_b__i_b_true = w_b__i_b_true;
    dataset.(testNames(k)).f_b__i_b_meas = f_b__i_b_meas;
    dataset.(testNames(k)).w_b__i_b_meas = w_b__i_b_meas;
    
end

