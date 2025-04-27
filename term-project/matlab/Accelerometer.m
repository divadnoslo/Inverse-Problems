classdef Accelerometer
    %ACCELEROMETER Summary of this class goes here
    %   Detailed explanation goes here
    

    %% Properties

    properties
        Bias (1,1) {isfloat, mustReal, mustBeFinite} = 0
        ScaleFactorError (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
        HigherOrderTerms (1,:) {isfloat, mustBeReal, mustBeFinite} = []
        AngleIntoPendulousAxis (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
        AngleIntoOutputAxis (1,1) {isfloat, mustBeReal, mustBeFinite} = 0
        
    end
    
    
    %% Methods

    % Constructor
    methods
        
        % Constructor
        function obj = Accelerometer()
        end
        
    end

    



end

