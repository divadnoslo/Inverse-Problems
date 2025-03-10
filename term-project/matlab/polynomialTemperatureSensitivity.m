function deviation = polynomialTemperatureSensitivity(coefficients, temperature, ambientTemperature)
% Polynomial Temperature Curve
%

%% Input and Output Arguments

arguments (Input)
    coefficients (1,:) {isfloat, mustBeReal, mustBeFinite}
    temperature (1,:) {isfloat, mustBeReal, mustBeFinite}
    ambientTemperature (1,1) {isfloat, mustBeReal, mustBeFinite}
end


%% Implementation

[~, K] = size(temperature);
deviation = coefficients(1) .* ones(1, K);

for k = 2 : length(coefficients)
    deviation = deviation + coefficients(k) * (temperature - ambientTemperature).^(k - 1);
end


end

