function C = Rz(psi)
% Rotate about Z
%
%

%% Input and Output Arguments

arguments (Input)
    psi (1,:) {isfloat, mustBeReal, mustBeFinite}
end

arguments (Output)
    C (3,3,:) {isfloat, mustBeReal, mustBeFinite}
end


%% Implementation

N = length(psi);
C = zeros(3, 3, N);

C(1,1,:) = cos(psi);
C(1,2,:) = -sin(psi);
C(2,1,:) = sin(psi);
C(2,2,:) = cos(psi);
C(3,3,:) = ones(1, N);

