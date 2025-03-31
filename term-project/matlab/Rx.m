function C = Rx(phi)
% Rotate about X
%
%

%% Input and Output Arguments

arguments (Input)
    phi (1,:) {isfloat, mustBeReal, mustBeFinite}
end

arguments (Output)
    C (3,3,:) {isfloat, mustBeReal, mustBeFinite}
end


%% Implementation

N = length(phi);
C = zeros(3, 3, N);

C(1,1,:) = ones(1, N);
C(2,2,:) = cos(phi);
C(2,3,:) = -sin(phi);
C(3,2,:) = sin(phi);
C(3,3,:) = cos(phi);

