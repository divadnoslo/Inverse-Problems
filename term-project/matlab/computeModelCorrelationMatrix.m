function Rho = computeModelCorrelationMatrix(C)
% Compute Model Correlation Matrix
%
%

%% Input and Output Arguments

arguments (Input)
    C (:,:) {isfloat, mustBeReal, mustBeFinite}
end

arguments (Output)
    Rho (:,:) {isfloat, mustBeReal, mustBeFinite}
end


%% Error Checking

[m, n] = size(C);
assert(...
    m == n, ...
    "computeModelCorrelationMatrix:matrixNotSquare", ...
    "Model covariance matrix must be square!")


%% Implemenation

Rho = zeros(n, n);
for ii = 1 : n
    for jj = 1 : n
        Rho(ii,jj) = C(ii,jj) ./ sqrt(C(ii,ii) * C(jj,jj));
    end
end


end

