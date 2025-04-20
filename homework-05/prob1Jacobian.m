function J = prob1Jacobian(m)
% Computes the Jacobian of eq. 9.56
%
% m = [A, f0, chi2, c]^T
%
%


%% Input and Output Arguments

arguments (Input)
    m (4,1) {isfloat, mustBeReal, mustBeFinite}
end

arguments (Output)
    J (2000,4) {isfloat, mustBeReal, mustBeFinite}
end


%% Implementation

persistent t
t = 0.02 * (0 : 1999).';

J = [...
    sin(2*pi*m(2).*t + m(3)), ...
    m(1)*2*pi.*t .* cos(2*pi*m(2).*t + m(3)), ...
    m(1) * cos(2*pi*m(2).*t + m(3)), ...
    ones(2000,1)];

