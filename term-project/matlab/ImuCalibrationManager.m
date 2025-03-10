classdef ImuCalibrationManager < handle
    % IMU Calibration Generator
    %
    %

    %% Properties

    properties (SetAccess = immutable)
        SamplingRate (1,1) {isfloat, mustBeReal, mustBePositive} = 100
        Duration (1,1) {isfloat, mustBeReal, mustBePositive} = 10
        Gravity (1,1) {isfloat, mustBeReal, mustBePositive} = 9.81
        SpinRate (1,1) {isfloat, mustBeReal, mustBeFinite} = 30 * pi/180
    end

    properties (Constant, Hidden)
        TestNames = [...
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
    end

    %% Methods

    % Constructor
    methods

        % Constructor
        function obj = ImuCalibrationManager(options)

            arguments (Input)
                options.SamplingRate (1,1) {isfloat, mustBeReal, mustBePositive} = 100
                options.Duration (1,1) {isfloat, mustBeReal, mustBePositive} = 10
                options.Gravity (1,1) {isfloat, mustBeReal, mustBePositive} = 9.81
                options.SpinRate (1,1) {isfloat, mustBeReal, mustBeFinite} = 30 * pi/180
            end

            obj.SamplingRate = options.SamplingRate;
            obj.Duration = options.SamplingRate;
            obj.Gravity = options.Gravity;
            obj.SpinRate = options.SpinRate;

        end

    end

    % Public Methods
    methods (Access = public)

        % Create Calibration Dataset
        function dataset = createCalibrationDataSet(obj, imu, temperature)

            arguments (Input)
                obj
                imu (1,1) ImuInterface
                temperature (1,1) {isfloat, mustBeReal, mustBeFinite} = 21
            end

            dataset = struct();

            % Simulate Tests
            dataset.time = 0 : (1/obj.SamplingRate) : obj.Duration;
            K = length(dataset.time);
            for k = 1 : length(obj.TestNames)

                % Initialize
                f_b__i_b_true = zeros(3, K);
                w_b__i_b_true = zeros(3, K);

                % Simulate Perfect Test
                switch obj.TestNames(k)

                    case "accelPosX"
                        f_b__i_b_true(1,:) = obj.Gravity * ones(1, K);

                    case "accelNegX"
                        f_b__i_b_true(1,:) = -1 * obj.Gravity * ones(1, K);

                    case "accelPosY"
                        f_b__i_b_true(2,:) = obj.Gravity * ones(1, K);

                    case "accelNegY"
                        f_b__i_b_true(2,:) = -1 * obj.Gravity * ones(1, K);

                    case "accelPosZ"
                        f_b__i_b_true(3,:) = obj.Gravity * ones(1, K);

                    case "accelNegZ"
                        f_b__i_b_true(3,:) = -1 * obj.Gravity * ones(1, K);

                    case "gyroPosX"
                        w_b__i_b_true(1,:) = obj.SpinRate * ones(1, K);

                    case "gyroNegX"
                        w_b__i_b_true(1,:) = -1 * obj.SpinRate * ones(1, K);

                    case "gyroPosY"
                        w_b__i_b_true(2,:) = obj.SpinRate * ones(1, K);

                    case "gyroNegY"
                        w_b__i_b_true(2,:) = -1 * obj.SpinRate * ones(1, K);

                    case "gyroPosZ"
                        w_b__i_b_true(3,:) = obj.SpinRate * ones(1, K);

                    case "gyroNegZ"
                        w_b__i_b_true(3,:) = -1 * obj.SpinRate * ones(1, K);

                    otherwise
                        error("Test not recognized!")

                end

                % Run IMU Forward Model
                [f_b__i_b_meas, w_b__i_b_meas] = imu.runForwardModel(...
                    f_b__i_b_true, ...
                    w_b__i_b_true);

                % Store Test Data
                dataset.(obj.TestNames(k)) = struct();
                dataset.(obj.TestNames(k)).f_b__i_b_true = f_b__i_b_true;
                dataset.(obj.TestNames(k)).w_b__i_b_true = w_b__i_b_true;
                dataset.(obj.TestNames(k)).f_b__i_b_meas = f_b__i_b_meas;
                dataset.(obj.TestNames(k)).w_b__i_b_meas = w_b__i_b_meas;

            end

        end

        % Plot Calibration Dataset
        function plotCalibrationDataset(obj, dataset)

            arguments (Input)
                obj
                dataset (1,1) {mustBeA(dataset, 'struct')}
            end

            XYZ = ["X", "Y", "Z"];

            % Accelerometer Dataset Plots
            for k = 1 : 3

                fig = figure("Name", XYZ(k) + "-Axis Up/Down Accelerometer Dataset");
                tl = tiledlayout(3, 1, "Parent", fig);

                ax = nexttile(1);
                hold(ax, "on")
                title("X-Axis Accelerometer Measurements")
                plot(dataset.time, dataset.("accelPos" + XYZ(k)).f_b__i_b_meas(1,:), 'r')
                plot(dataset.time, dataset.("accelNeg" + XYZ(k)).f_b__i_b_meas(1,:), 'k')
                xlabel("Time [sec]")
                ylabel("[m/sec/sec]")
                grid on
                grid minor

                ax = nexttile(2);
                hold(ax, "on")
                title("Y-Axis Accelerometer Measurements")
                plot(dataset.time, dataset.("accelPos" + XYZ(k)).f_b__i_b_meas(2,:), 'g')
                plot(dataset.time, dataset.("accelNeg" + XYZ(k)).f_b__i_b_meas(2,:), 'k')
                xlabel("Time [sec]")
                ylabel("[m/sec/sec]")
                grid on
                grid minor

                ax = nexttile(3);
                hold(ax, "on")
                title("Z-Axis Accelerometer Measurements")
                plot(dataset.time, dataset.("accelPos" + XYZ(k)).f_b__i_b_meas(3,:), 'b')
                plot(dataset.time, dataset.("accelNeg" + XYZ(k)).f_b__i_b_meas(3,:), 'k')
                xlabel("Time [sec]")
                ylabel("[m/sec/sec]")
                grid on
                grid minor

                linkaxes(tl.Children, 'x')

            end

            % Gyroscope Dataset Plots
            for k = 1 : 3

                fig = figure("Name", XYZ(k) + "-Axis Pos/Neg Gyroscope Dataset");
                tl = tiledlayout(3, 1, "Parent", fig);

                ax = nexttile(1);
                hold(ax, "on")
                title("X-Axis Gyroscope Measurements")
                plot(dataset.time, 180/pi * dataset.("gyroPos" + XYZ(k)).w_b__i_b_meas(1,:), 'r')
                plot(dataset.time, 180/pi * dataset.("gyroNeg" + XYZ(k)).w_b__i_b_meas(1,:), 'k')
                xlabel("Time [sec]")
                ylabel("[deg/sec]")
                grid on
                grid minor

                ax = nexttile(2);
                hold(ax, "on")
                title("Y-Axis Gyroscope Measurements")
                plot(dataset.time, 180/pi * dataset.("gyroPos" + XYZ(k)).w_b__i_b_meas(2,:), 'g')
                plot(dataset.time, 180/pi * dataset.("gyroNeg" + XYZ(k)).w_b__i_b_meas(2,:), 'k')
                xlabel("Time [sec]")
                ylabel("[deg/sec]")
                grid on
                grid minor

                ax = nexttile(3);
                hold(ax, "on")
                title("Z-Axis Gyroscope Measurements")
                plot(dataset.time, 180/pi * dataset.("gyroPos" + XYZ(k)).w_b__i_b_meas(3,:), 'b')
                plot(dataset.time, 180/pi * dataset.("gyroNeg" + XYZ(k)).w_b__i_b_meas(3,:), 'k')
                xlabel("Time [sec]")
                ylabel("[degsec]")
                grid on
                grid minor

                linkaxes(tl.Children, 'x')

            end

        end
        
        function imu = processCalibrationDataSet(obj, dataset)

            arguments (Input)
                obj
                dataset (1,1) {mustBeA(dataset, 'struct')}
            end

            XYZ = ["X", "Y", "Z"];
            XYZM = ["XY", "XZ", "YX", "YZ", "ZX", "ZY"];

            imu = ImuModel();

            % Accelerometer Fixed Bias
            for k = 1 : 3
                imu.("AccelFixedBias" + XYZ(k)) = ...
                    (1/2) * mean(...
                    dataset.("accelPos" + XYZ(k)).f_b__i_b_meas(k,:) + dataset.("accelNeg" + XYZ(k)).f_b__i_b_meas(k,:));
            end

            % Accelerometer Scale Factor Error
            for k = 1 : 3
                imu.("AccelScaleFactorError" + XYZ(k)) = ...
                    (1 / (2 * obj.Gravity)) * mean(...
                    dataset.("accelPos" + XYZ(k)).f_b__i_b_meas(k,:) - dataset.("accelNeg" + XYZ(k)).f_b__i_b_meas(k,:)) ...
                    - 1;
            end

            % Accelerometer Misalignment
            for k = 1 : 6
                
                switch XYZM(k)
                    
                    case "XY"
                        imu.("AccelMisalignment" + XYZM(k)) = ...
                            (1 / (2 * obj.Gravity)) * mean(...
                            dataset.("accelPos" + XYZ(2)).f_b__i_b_meas(1,:) - dataset.("accelNeg" + XYZ(2)).f_b__i_b_meas(1,:));

                    case "XZ"
                        imu.("AccelMisalignment" + XYZM(k)) = ...
                            (1 / (2 * obj.Gravity)) * mean(...
                            dataset.("accelPos" + XYZ(3)).f_b__i_b_meas(1,:) - dataset.("accelNeg" + XYZ(3)).f_b__i_b_meas(1,:));

                    case "YX"
                        imu.("AccelMisalignment" + XYZM(k)) = ...
                            (1 / (2 * obj.Gravity)) * mean(...
                            dataset.("accelPos" + XYZ(1)).f_b__i_b_meas(2,:) - dataset.("accelNeg" + XYZ(1)).f_b__i_b_meas(2,:));

                    case "YZ"
                        imu.("AccelMisalignment" + XYZM(k)) = ...
                            (1 / (2 * obj.Gravity)) * mean(...
                            dataset.("accelPos" + XYZ(3)).f_b__i_b_meas(2,:) - dataset.("accelNeg" + XYZ(3)).f_b__i_b_meas(2,:));

                    case "ZX"
                        imu.("AccelMisalignment" + XYZM(k)) = ...
                            (1 / (2 * obj.Gravity)) * mean(...
                            dataset.("accelPos" + XYZ(1)).f_b__i_b_meas(3,:) - dataset.("accelNeg" + XYZ(1)).f_b__i_b_meas(3,:));

                    case "ZY"
                        imu.("AccelMisalignment" + XYZM(k)) = ...
                            (1 / (2 * obj.Gravity)) * mean(...
                            dataset.("accelPos" + XYZ(2)).f_b__i_b_meas(3,:) - dataset.("accelNeg" + XYZ(2)).f_b__i_b_meas(3,:));

                    otherwise
                        error("Misalignment case not recognized!")

                end

            end

            % Gyroscope Fixed Bias
            for k = 1 : 3
                imu.("GyroFixedBias" + XYZ(k)) = ...
                    (1/2) * mean(...
                    dataset.("gyroPos" + XYZ(k)).w_b__i_b_meas(k,:) + dataset.("gyroNeg" + XYZ(k)).w_b__i_b_meas(k,:));
            end

            % Gyroscope Scale Factor Error
            for k = 1 : 3
                imu.("GyroScaleFactorError" + XYZ(k)) = ...
                    (1 / (2 * obj.SpinRate)) * mean(...
                    dataset.("gyroPos" + XYZ(k)).w_b__i_b_meas(k,:) - dataset.("gyroNeg" + XYZ(k)).w_b__i_b_meas(k,:)) ...
                    - 1;
            end

            % Gyroscope Misalignment
            for k = 1 : 6
                
                switch XYZM(k)
                    
                    case "XY"
                        imu.("GyroMisalignment" + XYZM(k)) = ...
                            (1 / (2 * obj.SpinRate)) * mean(...
                            dataset.("gyroPos" + XYZ(2)).w_b__i_b_meas(1,:) - dataset.("gyroNeg" + XYZ(2)).w_b__i_b_meas(1,:));

                    case "XZ"
                        imu.("GyroMisalignment" + XYZM(k)) = ...
                            (1 / (2 * obj.SpinRate)) * mean(...
                            dataset.("gyroPos" + XYZ(3)).w_b__i_b_meas(1,:) - dataset.("gyroNeg" + XYZ(3)).w_b__i_b_meas(1,:));

                    case "YX"
                        imu.("GyroMisalignment" + XYZM(k)) = ...
                            (1 / (2 * obj.SpinRate)) * mean(...
                            dataset.("gyroPos" + XYZ(1)).w_b__i_b_meas(2,:) - dataset.("gyroNeg" + XYZ(1)).w_b__i_b_meas(2,:));

                    case "YZ"
                        imu.("GyroMisalignment" + XYZM(k)) = ...
                            (1 / (2 * obj.SpinRate)) * mean(...
                            dataset.("gyroPos" + XYZ(3)).w_b__i_b_meas(2,:) - dataset.("gyroNeg" + XYZ(3)).w_b__i_b_meas(2,:));

                    case "ZX"
                        imu.("GyroMisalignment" + XYZM(k)) = ...
                            (1 / (2 * obj.SpinRate)) * mean(...
                            dataset.("gyroPos" + XYZ(1)).w_b__i_b_meas(3,:) - dataset.("gyroNeg" + XYZ(1)).w_b__i_b_meas(3,:));

                    case "ZY"
                        imu.("GyroMisalignment" + XYZM(k)) = ...
                            (1 / (2 * obj.SpinRate)) * mean(...
                            dataset.("gyroPos" + XYZ(2)).w_b__i_b_meas(3,:) - dataset.("gyroNeg" + XYZ(2)).w_b__i_b_meas(3,:));

                    otherwise
                        error("Misalignment case not recognized!")

                end

            end

        end

    end

end

