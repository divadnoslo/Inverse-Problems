classdef NormalDistribution
    % Normal Distribution
    
    %% Properties

    properties
        Mean (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
        Variance (1,1) {isfloat, mustBeReal, mustBeFinite, mustBePositive} = 1
    end

    properties (Dependent)
        StandardDeviation (1,1) {isfloat, mustBeReal, mustBeFinite, mustBePositive}
    end
    
    %% Methods

    % Constructor, Set/Get Methods
    methods

        % Constructor
        function obj = NormalDistribution(mean, variance)

            arguments (Input)
                mean (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
                variance (1,1) {isfloat, mustBeReal, mustBeFinite, mustBePositive} = 1
            end

            obj.Mean = mean;
            obj.Variance = variance;

        end
        
        % Set Standard Deviation
        function obj = set.StandardDeviation(obj, standardDeviation)

            arguments (Input)
                obj
                standardDeviation (1,1) {isfloat, mustBeReal, mustBeFinite, mustBePositive}
            end

            obj.Variance = standardDeviation^2;

        end

        % Get Standard Deviation
        function standardDeviation = get.StandardDeviation(obj)
            standardDeviation = sqrt(obj.Variance);
        end

    end

    % Public Methods
    methods (Access = public)

        function p = probabilityDensityFunction(obj, x)

            arguments (Input)
                obj
                x (1,:) {isfloat, mustBeReal, mustBeFinite}
            end

            arguments (Output)
                p (1,:) {isfloat, mustBeReal, mustBeFinite}
            end

            p = ...
                (1 / (sqrt(2*pi) * obj.StandardDeviation)) * ...
                exp((-1/2) * (x - obj.Mean).^2 ./ obj.Variance);

        end

    end

end

