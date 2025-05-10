function L = chi2likelihood(chi2)
% Compute Likelihood Given a Chi^2 Value
%
%

%% Input and Output Arguments

arguments (Input)
    chi2 (:,:) {isfloat, mustBeReal, mustBeFinite}
end

arguments (Output)
    L (:,:) {isfloat, mustBeReal, mustBeFinite}
end

%% Implementation

% Persistent Call
persistent be10 const
be10 = load(fullfile(pwd, "data", "be10.mat"));
const = prod(1 ./ (be10.sigma .* sqrt(2*pi)));

% Compute Likelihood
L = const * exp(-chi2/2);

end

