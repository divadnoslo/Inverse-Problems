classdef Wgs84
    % A class to hold WGS-84 Parameters
    
    properties (Constant)
        R0 = 6378137;
        Rp = 6356752.31425;
        f = 1 / 298.257223563;
        e = 0.0818191908425;
    end
    
    methods

        function obj = Wgs84()
        end

    end

end

