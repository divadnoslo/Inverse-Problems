% Prob2 Uniform Log-Prior
%
function lp = logprior2(m)

persistent N
N = NormalDistribution(0.0005, 0.0002^2);
lp = log(N.probabilityDensityFunction(m(1)));


