function N = func(m)
% Compute Concentration (N)
%
% m = [ep, T]
%
%

%% Input and Output Arguments

arguments (Input)
    m (2,1) {isfloat, mustBeReal, mustBeFinite}
end

arguments (Output)
    N (:,1) {isfloat, mustBeReal}
end


%% Implementation

persistent lambda P0 mu d
lambda = 4.998e-7;
P0 = 3.2;
mu = 0.0166;
d = load(fullfile(pwd, "data", "be10.mat"));

N = (P0/(lambda + mu*m(1))) .* exp(-mu.*d.depths) .* (1 - exp(-1 * (lambda + mu*m(1)) .* m(2)));


end
