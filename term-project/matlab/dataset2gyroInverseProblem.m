function [G, d] = dataset2gyroInverseProblem(dataset)
% Create a Model Operator from a Traditional Calibration Dataset
%
%


%% Input and Output Arguments

arguments (Input)
    dataset (1,1) struct
end

%% Implementation

K = length(dataset.time);

% Initialize
m = 3 * K;
n = 4 * 3;
G = zeros(m, n);
d = zeros(m, 1);

% Build from X Tests
G((1:K) + 0*K, (1:4) + 0*4) = [ones(K, 1), dataset.gyroPosX.w_b__i_b_true.'];
G((1:K) + 1*K, (1:4) + 0*4) = [ones(K, 1), dataset.gyroNegX.w_b__i_b_true.'];
d((1:K) + 0*K) = dataset.gyroPosX.w_b__i_b_meas(1,:).' - dataset.gyroPosX.w_b__i_b_true(1,:).';
d((1:K) + 1*K) = dataset.gyroNegX.w_b__i_b_meas(1,:).' - dataset.gyroNegX.w_b__i_b_true(1,:).';

% Build from Y Tests
G((1:K) + 2*K, (1:4) + 1*4) = [ones(K, 1), dataset.gyroPosY.w_b__i_b_true.'];
G((1:K) + 3*K, (1:4) + 1*4) = [ones(K, 1), dataset.gyroNegY.w_b__i_b_true.'];
d((1:K) + 2*K) = dataset.gyroPosY.w_b__i_b_meas(2,:).' - dataset.gyroPosY.w_b__i_b_true(2,:).';
d((1:K) + 3*K) = dataset.gyroNegY.w_b__i_b_meas(2,:).' - dataset.gyroNegY.w_b__i_b_true(2,:).';

% Build from Z Tests
G((1:K) + 4*K, (1:4) + 2*4) = [ones(K, 1), dataset.gyroPosZ.w_b__i_b_true.'];
G((1:K) + 5*K, (1:4) + 2*4) = [ones(K, 1), dataset.gyroNegZ.w_b__i_b_true.'];
d((1:K) + 4*K) = dataset.gyroPosZ.w_b__i_b_meas(3,:).' - dataset.gyroPosZ.w_b__i_b_true(3,:).';
d((1:K) + 5*K) = dataset.gyroNegZ.w_b__i_b_meas(3,:).' - dataset.gyroNegZ.w_b__i_b_true(3,:).';


end

