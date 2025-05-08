function f = log_of_f_d_given_m(m)
% Compute Likelihood Given a Model
%
%

%% Input and Output Arguments

arguments (Input)
    m (2,1) {isfloat, mustBeReal, mustBeFinite}
end

arguments (Output)
    f (1,1) % {isfloat, mustBeReal, mustBeFinite}
end

%% Implementation

% Persistent Call
persistent be10 K a b
be10 = load(fullfile(pwd, "data", "be10.mat"));
K = length(be10.d);
a = 1 ./ (be10.sigma*sqrt(2*pi));
b = 1 ./ (2 * be10.sigma.^2);

% Compute Likelihood
f = sum((func(m) - be10.d).^2 ./ (2 * be10.sigma.^2));

dummy = 1;

end

