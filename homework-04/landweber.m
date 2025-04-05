function X = landweber(G, m0, d, omega, k)
% Landweber Iteration for Sttepest Decent Method
%
% Performs "k" iterations of the Landweber interation on Gm = d, with
% initial solution m0, and parameter omega.
%
% On output, X is a matrix of size n by k. The k columns of X are the
% iterates m1, m2, ..., mk.
%
% For convergence, we need 0 < omega < 2/s(1)^2. Outside of this range, you
% can expect the method to diverge.

%% Input and Output Arguments

arguments (Input)
    G (:,:) {isfloat, mustBeReal, mustBeFinite}
    m0 (:,1) {isfloat, mustBeReal, mustBeFinite}
    d (:,1) {isfloat, mustBeReal, mustBeFinite}
    omega (1,1) {isfloat, mustBeReal, mustBeFinite, mustBePositive}
    k (1,1) {mustBeInteger, mustBePositive}
end

arguments (Output)
    X (:,:) {isfloat, mustBeReal, mustBeFinite}
end


%% Error Checking

[m, n] = size(G);
assert(...
    n == length(m0), ...
    "Initial model m0 must have the same number of columns as G!");
assert(...
    m == length(d), ...
    "Data vector d must have the same number of rows as G!");

assert(...
    omega < (2/(svds(G,1)^2)), ...
    "This omega value will not let the iterations converge to a solution!")


%% Implementation

X = zeros(n, k);
X(:,1) = m0;
for ii = 2 : k
    X(:,ii) = X(:,ii-1) - omega * G.' * (G*X(:,ii-1) - d);
end


end

