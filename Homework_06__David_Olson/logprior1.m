% Prob2 Uniform Log-Prior
%
function lp = logprior1(m)

if ((m(1) >= 5e-6) && (m(1) <= 1e-3)) && ((m(2) >= 500) && (m(2) <= 199500))
  lp = 0;
else
  lp = -Inf;
end
