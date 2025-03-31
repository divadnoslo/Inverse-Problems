function C = Ry(theta)
% Rotate about Y
%
%

%% Input and Output Arguments

arguments (Input)
    theta (1,:) {isfloat, mustBeReal, mustBeFinite}
end

arguments (Output)
    C (3,3,:) {isfloat, mustBeReal, mustBeFinite}
end


%% Implementation

N = length(theta);
C = zeros(3, 3, N);

C(1,1,:) = cos(theta);
C(1,3,:) = sin(theta);
C(2,2,:) = ones(1, N);
C(3,1,:) = -sin(theta);
C(3,3,:) = cos(theta);

