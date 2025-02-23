function fillHandle = plotConfidenceEllipse(modelParameters, modelCovariance, degreesOfFreedom, confidenceInterval, options)
% Plot Confidence Ellipse
%
% This function plots a 2D ellipse based upon a set of model parameters,
% model covariance, degrees of freedom in the Chi-Squared probability
% density function, and the confidence interval region of interest.
%
% Inputs
% - modelParameters
%     - The weighted least squares solution
% - modelCovariance
%     - The model covariance matrix
% - degreesOfFreedom
%     - The degrees of freedom as defined in the Chi-Squared distribution
% - confidenceInterval
%     - The confidence percentage expressed in decimal form
%     - Default: 95%
%
% Optional Name-Value Pair Arguments
% - NumberOfPoints
%     - The number of points used to plot the ellipse
%     - Default: 100
% - Color
%     - The color of the plotted ellipse
%     - Default: blue
% - FaceAlpha
%     - The transparency of the ellipse
%     - Default: 0.4
% - Axes
%     - The axes of the ellipse to be plotted upon
%     - Default: gca
%
% References:
% Parameter Estimation and Inverse Problems, 3rd Edition
% By Richard C. Aster, Brian Borchers, and Clifford H. Thurber
%
% https://github.com/brianborchers/PEIP/tree/master/Examples/chap2/ex_2_1_2


%% Input and Output Arguments

arguments (Input)
    modelParameters (2,1) {isfloat, mustBeReal, mustBeFinite}
    modelCovariance (2,2) {isfloat, mustBeReal, mustBeFinite}
    degreesOfFreedom (1,1) {isfloat, mustBeReal, mustBeFinite, mustBePositive, mustBeInteger}
    confidenceInterval (1,1) {isfloat, mustBeReal, mustBeInRange(confidenceInterval, 0, 1, "exclusive")} = 0.95
    options.NumberOfPoints (1,1) {mustBeReal, mustBeFinite, mustBePositive} = 100
    options.Color (1,1) string {mustBeText} = "b"
    options.FaceAlpha (1,1) {isfloat, mustBeReal, mustBeInRange(options.FaceAlpha, 0, 1)} = 0.4
    options.Axes matlab.graphics.axis.Axes = gca
end

arguments (Output)
    fillHandle matlab.graphics.primitive.Patch
end


%% Error Checking

% Assert model covariance is a symmetric matrix
assert(...
    all(modelCovariance == modelCovariance.', 'all'), ...
    "peip:plotConfidenceEllipse:modelCovarianceIsNotSymmetric", ...
    "Provided model covariance matrix must be symmetric!")

% Assert model covariance is positive definite
assert(...
    trace(modelCovariance) > 0, ...
    "peip:plotConfidenceEllipse:modelCovarianceTraceIsNotPositive", ...
    "To be positive semi-definite, the trace of the model covariance must be positive!")


%% Implementation

% Compute confidence interval region
delta2 = chi2inv(confidenceInterval, degreesOfFreedom);
delta = sqrt(delta2);

% Compute ellipse semi-major and semi-minor axes
[U, Lambda] = eig(inv(modelCovariance));

% Compute points for ellipse plotting
theta = linspace(0, 2*pi, options.NumberOfPoints);
radius = zeros(options.NumberOfPoints, 2);
radius(:,1) = ...
    (delta / sqrt(Lambda(1,1))) * U(1,1) * cos(theta) + ...
    (delta / sqrt(Lambda(2,2))) * U(1,2) * sin(theta);
radius(:,2) = ...
    (delta / sqrt(Lambda(1,1))) * U(2,1) * cos(theta) + ...
    (delta / sqrt(Lambda(2,2))) * U(2,2) * sin(theta);

% Plot ellipse
fillHandle = fill(...
    options.Axes, ...
    modelParameters(1) + radius(:,1), ...
    modelParameters(2) + radius(:,2), ...
    options.Color, ...
    "FaceAlpha", options.FaceAlpha);


end

