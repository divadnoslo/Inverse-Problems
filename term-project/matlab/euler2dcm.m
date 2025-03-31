function C = euler2dcm(phi, theta, psi)
% Convert Euler Angle Rotation Sequence to DCM
%
%

%% Input and Output Arguments

arguments (Input)
    phi (1,:) {isfloat, mustBeReal, mustBeFinite}
    theta (1,:) {isfloat, mustBeReal, mustBeFinite}
    psi (1,:) {isfloat, mustBeReal, mustBeFinite}
end

arguments (Output)
    C (3,3,:) {isfloat, mustBeReal, mustBeFinite}
end


%% Error Checking

N = length(phi);
assert(...
    N == length(theta), ...
    "The length of theta must match the length of phi!")
assert(...
    N == length(psi), ...
    "The length of psi must match the length of phi!")


%% Implementation

C = pagemtimes(Rz(psi), pagemtimes(Ry(theta), Rx(phi)));


end

