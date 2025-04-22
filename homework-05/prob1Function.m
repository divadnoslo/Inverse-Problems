function F = prob1Function(m)
% Implements the function in eq. 9.56
%
% m = [A, f0, chi2, sc^T
%
%

%% Input and Output Arguments

arguments (Input)
    m (:,1) {isfloat, mustBeReal, mustBeFinite}
end

arguments (Output)
    F (:,1) {isfloat, mustBeReal, mustBeFinite}
end


%% Implementation

persistent t
t = 0.02 * (0 : 1999).';

persistent d
d = load(fullfile(pwd, "data", "instdata.mat"));

F = (m(1)*sin(2*pi*m(2).*t + m(3)) + m(4)) - d.instdata;


end

