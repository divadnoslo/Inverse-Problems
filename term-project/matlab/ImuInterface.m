classdef (Abstract) ImuInterface < handle
    % IMU Interface Abstract Class to define the necessary properties and
    % methods that every IMU related calss should have
    
    methods (Abstract)
        runForwardModel(obj)
        runInverseModel(obj)
    end
end

