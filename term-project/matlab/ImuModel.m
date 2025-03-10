classdef ImuModel < ImuInterface
    % IMU Model as decribed by Groves
    
    %% Properties

    % Accelerometer Properties
    properties
        AccelFixedBiasX (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
        AccelFixedBiasY (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
        AccelFixedBiasZ (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
        AccelScaleFactorErrorX (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
        AccelScaleFactorErrorY (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
        AccelScaleFactorErrorZ (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
        AccelMisalignmentXY (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
        AccelMisalignmentXZ (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
        AccelMisalignmentYX (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
        AccelMisalignmentYZ (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
        AccelMisalignmentZX (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
        AccelMisalignmentZY (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
        AccelWhiteNoiseOneSigmaX (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
        AccelWhiteNoiseOneSigmaY (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
        AccelWhiteNoiseOneSigmaZ (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
    end
    
    % Accelerometer Dependent Properties
    properties (Dependent, Hidden)
        b_a (3,1) {isfloat, mustBeReal, mustBeFinite}
        M_a (3,3) {isfloat, mustBeReal, mustBeFinite}
        sigma_a (3,1) {isfloat, mustBeReal, mustBeFinite}
    end

    % Gyroscope Properties
    properties
        GyroFixedBiasX (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
        GyroFixedBiasY (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
        GyroFixedBiasZ (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
        GyroScaleFactorErrorX (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
        GyroScaleFactorErrorY (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
        GyroScaleFactorErrorZ (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
        GyroMisalignmentXY (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
        GyroMisalignmentXZ (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
        GyroMisalignmentYX (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
        GyroMisalignmentYZ (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
        GyroMisalignmentZX (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
        GyroMisalignmentZY (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
        GyroWhiteNoiseOneSigmaX (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
        GyroWhiteNoiseOneSigmaY (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
        GyroWhiteNoiseOneSigmaZ (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
    end
    
    % Gyroscope Dependent Properties
    properties (Dependent, Hidden)
        b_g (3,1) {isfloat, mustBeReal, mustBeFinite}
        M_g (3,3) {isfloat, mustBeReal, mustBeFinite}
        sigma_g (3,1) {isfloat, mustBeReal, mustBeFinite}
    end

    %% Methods

    % Constructor
    methods

        % Constructor
        function obj = BasicImuModel()
            % Basic IMU Model
        end
        
    end

    % Accelerometer Set/Get Methods
    methods

        % Set Accelerometer Fixed Bias
        function set.b_a(obj, b_a)

            arguments (Input)
                obj
                b_a (3,1) {isfloat, mustBeReal, mustBeFinite}
            end

            obj.AccelFixedBiasX = b_a(1);
            obj.AccelFixedBiasY = b_a(2);
            obj.AccelFixedBiasZ = b_a(3);

        end

        % Get Accelerometer Fixed Bias
        function b = get.b_a(obj)

            arguments (Output)
                b (3,1) {isfloat, mustBeReal, mustBeFinite}
            end

            b = zeros(3, 1);
            b(1) = obj.AccelFixedBiasX;
            b(2) = obj.AccelFixedBiasY;
            b(3) = obj.AccelFixedBiasZ;

        end

        % Set Accelerometer Scale Factor and Misalignment
        function set.M_a(obj, M)

            arguments (Input)
                obj
                M (3,3) {isfloat, mustBeReal, mustBeFinite}
            end

            obj.AccelScaleFactorErrorX = M(1,1);
            obj.AccelMisalignmentXY = M(1,2);
            obj.AccelMisalignmentXZ = M(1,3);

            obj.AccelMisalignmentYX = M(2,1);
            obj.AccelScaleFactorErrorY = M(2,2);
            obj.AccelMisalignmentYZ = M(2,3);

            obj.AccelMisalignmentZX = M(3,1);
            obj.AccelMisalignmentZY = M(3,2);
            obj.AccelScaleFactorErrorZ = M(3,3);

        end

        % Get Accelerometer Scale Factor and Misalignment
        function M = get.M_a(obj)

            arguments (Output)
                M (3,3) {isfloat, mustBeReal, mustBeFinite}
            end

            M = zeros(3, 3);

            M(1,1) = obj.AccelScaleFactorErrorX;
            M(1,2) = obj.AccelMisalignmentXY;
            M(1,3) = obj.AccelMisalignmentXZ;

            M(2,1) = obj.AccelMisalignmentYX;
            M(2,2) = obj.AccelScaleFactorErrorY;
            M(2,3) = obj.AccelMisalignmentYZ;

            M(3,1) = obj.AccelMisalignmentZX;
            M(3,2) = obj.AccelMisalignmentZY;
            M(3,3) = obj.AccelScaleFactorErrorZ;

        end

        % Set Accelerometer "sigma_a"
        function set.sigma_a(obj, sigmaVector)

            arguments (Input)
                obj
                sigmaVector (3,1) {isfloat, mustBeReal, mustBeFinite}
            end

            obj.AccelWhiteNoiseOneSigmaX = sigmaVector(1);
            obj.AccelWhiteNoiseOneSigmaY = sigmaVector(2);
            obj.AccelWhiteNoiseOneSigmaZ = sigmaVector(3);

        end

        % Get Accelerometer "sigma_a"
        function sigmaVector = get.sigma_a(obj)

            arguments (Output)
                sigmaVector (3,1) {isfloat, mustBeReal, mustBeFinite}
            end

            sigmaVector = zeros(3, 1);
            sigmaVector(1) = obj.AccelWhiteNoiseOneSigmaX;
            sigmaVector(2) = obj.AccelWhiteNoiseOneSigmaY;
            sigmaVector(3) = obj.AccelWhiteNoiseOneSigmaZ;

        end

    end

    % Gyroscope Set/Get Methods
    methods

        % Set Gyroscope Fixed Bias
        function set.b_g(obj, b_g)

            arguments (Input)
                obj
                b_g (3,1) {isfloat, mustBeReal, mustBeFinite}
            end

            obj.GyroFixedBiasX = b_g(1);
            obj.GyroFixedBiasY = b_g(2);
            obj.GyroFixedBiasZ = b_g(3);

        end

        % Get Gyroscope Fixed Bias
        function b = get.b_g(obj)

            arguments (Output)
                b (3,1) {isfloat, mustBeReal, mustBeFinite}
            end

            b = zeros(3, 1);
            b(1) = obj.GyroFixedBiasX;
            b(2) = obj.GyroFixedBiasY;
            b(3) = obj.GyroFixedBiasZ;

        end

        % Set Gyroscope Scale Factor and Misalignment
        function set.M_g(obj, M)

            arguments (Input)
                obj
                M (3,3) {isfloat, mustBeReal, mustBeFinite}
            end

            obj.GyroScaleFactorErrorX = M(1,1);
            obj.GyroMisalignmentXY = M(1,2);
            obj.GyroMisalignmentXZ = M(1,3);

            obj.GyroMisalignmentYX = M(2,1);
            obj.GyroScaleFactorErrorY = M(2,2);
            obj.GyroMisalignmentYZ = M(2,3);

            obj.GyroMisalignmentZX = M(3,1);
            obj.GyroMisalignmentZY = M(3,2);
            obj.GyroScaleFactorErrorZ = M(3,3);

        end

        % Get Gyroscope Scale Factor and Misalignment
        function M = get.M_g(obj)

            arguments (Output)
                M (3,3) {isfloat, mustBeReal, mustBeFinite}
            end

            M = zeros(3, 3);

            M(1,1) = obj.GyroScaleFactorErrorX;
            M(1,2) = obj.GyroMisalignmentXY;
            M(1,3) = obj.GyroMisalignmentXZ;

            M(2,1) = obj.GyroMisalignmentYX;
            M(2,2) = obj.GyroScaleFactorErrorY;
            M(2,3) = obj.GyroMisalignmentYZ;

            M(3,1) = obj.GyroMisalignmentZX;
            M(3,2) = obj.GyroMisalignmentZY;
            M(3,3) = obj.GyroScaleFactorErrorZ;

        end

        % Set Gyroscope "sigma_g"
        function obj = set.sigma_g(obj, sigmaVector)

            arguments (Input)
                obj
                sigmaVector (3,1) {isfloat, mustBeReal, mustBeFinite}
            end

            obj.GyroWhiteNoiseOneSigmaX = sigmaVector(1);
            obj.GyroWhiteNoiseOneSigmaY = sigmaVector(2);
            obj.GyroWhiteNoiseOneSigmaZ = sigmaVector(3);

        end

        % Get Gyroscope "sigma_g"
        function sigmaVector = get.sigma_g(obj)

            arguments (Output)
                sigmaVector (3,1) {isfloat, mustBeReal, mustBeFinite}
            end

            sigmaVector = zeros(3, 1);
            sigmaVector(1) = obj.GyroWhiteNoiseOneSigmaX;
            sigmaVector(2) = obj.GyroWhiteNoiseOneSigmaY;
            sigmaVector(3) = obj.GyroWhiteNoiseOneSigmaZ;

        end

    end

    % Public Methods
    methods (Access = public)

        % Run Forward Model
        function [f_b__i_b_meas, w_b__i_b_meas] = runForwardModel(obj, f_b__i_b_true, w_b__i_b_true)

            arguments (Input)
                obj
                f_b__i_b_true (3,:) {isfloat, mustBeReal, mustBeFinite}
                w_b__i_b_true (3,:) {isfloat, mustBeReal, mustBeFinite}
            end

            arguments (Output)
                f_b__i_b_meas (3,:) {isfloat, mustBeReal, mustBeFinite}
                w_b__i_b_meas (3,:) {isfloat, mustBeReal, mustBeFinite}
            end

            [~, K1] = size(f_b__i_b_true);
            [~, K2] = size(w_b__i_b_true);

            assert(...
                K1 == K2, ...
                "basicImuModel:runForwardModel:dataLengthMismatch", ...
                "Specific force and angular velocity truth must have the same number of columns!")

            f_b__i_b_meas = (eye(3) + obj.M_a) * f_b__i_b_true + obj.b_a + ...
                obj.sigma_a .* randn(3, K1);

            w_b__i_b_meas = (eye(3) + obj.M_g) * w_b__i_b_true + obj.b_g + ...
                obj.sigma_g .* randn(3, K2);

        end

        % Run Inverse Model
        function [f_b__i_b_comp, w_b__i_b_comp] = runInverseModel(obj, f_b__i_b_meas, w_b__i_b_meas)

            arguments (Input)
                obj
                f_b__i_b_meas (3,:) {isfloat, mustBeReal, mustBeFinite}
                w_b__i_b_meas (3,:) {isfloat, mustBeReal, mustBeFinite}
            end

            arguments (Output)
                f_b__i_b_comp (3,:) {isfloat, mustBeReal, mustBeFinite}
                w_b__i_b_comp (3,:) {isfloat, mustBeReal, mustBeFinite}
            end

            [~, K1] = size(f_b__i_b_meas);
            [~, K2] = size(w_b__i_b_meas);

            assert(...
                K1 == K2, ...
                "basicImuModel:runInverseModel:dataLengthMismatch", ...
                "Specific force and angular velocity measurements must have the same number of columns!")

            f_b__i_b_comp = (eye(3) + obj.M_a) \ (f_b__i_b_meas - obj.b_a);

            w_b__i_b_comp = (eye(3) + obj.M_g) \ (w_b__i_b_meas - obj.b_g);

        end

    end

end

