function y = prob1Function(m)
% Implements the function in eq. 9.56
%
% m = [A, f0, c, s]^T
%
%

%% Input and Output Arguments

arguments (Input)
    m (:,1) {isfloat, mustBeReal, mustBeFinite}
end

arguments (Output)
    y (:,1) {isfloat, mustBeReal, mustBeFinite}
end


%% Implementation

persistent t
t = 0.02 * (0 : 1999).';

y = m(1)*sin(2*pi*m(2).*t + m(3)) + m(4);


end

