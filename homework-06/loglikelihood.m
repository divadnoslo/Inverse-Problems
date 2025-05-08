% Example 11.4
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% l=loglikelihood(m)
%
function l = loglikelihood(m)

persistent lambda P0 mu d
lambda = 4.998e-7;
P0 = 3.2;
mu = 0.0166;
d = load(fullfile(pwd, "data", "be10.mat"));

N = (P0/(lambda + mu*m(1))) .* exp(-mu.*d.depths) .* (1 - exp(-1 * (lambda + mu*m(1)) .* m(2)));
chi2 = sum((N - d.d).^2 ./ (d.sigma.^2)); 
l = (-1/2) * chi2;