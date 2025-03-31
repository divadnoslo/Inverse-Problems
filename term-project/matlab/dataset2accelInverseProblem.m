function [G, d] = dataset2accelInverseProblem(dataset)
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
G((1:K) + 0*K, (1:4) + 0*4) = [ones(K, 1), dataset.accelPosX.f_b__i_b_true.'];
G((1:K) + 1*K, (1:4) + 0*4) = [ones(K, 1), dataset.accelNegX.f_b__i_b_true.'];
d((1:K) + 0*K) = dataset.accelPosX.f_b__i_b_meas(1,:).' - dataset.accelPosX.f_b__i_b_true(1,:).';
d((1:K) + 1*K) = dataset.accelNegX.f_b__i_b_meas(1,:).' - dataset.accelNegX.f_b__i_b_true(1,:).';

% Build from Y Tests
G((1:K) + 2*K, (1:4) + 1*4) = [ones(K, 1), dataset.accelPosY.f_b__i_b_true.'];
G((1:K) + 3*K, (1:4) + 1*4) = [ones(K, 1), dataset.accelNegY.f_b__i_b_true.'];
d((1:K) + 2*K) = dataset.accelPosY.f_b__i_b_meas(2,:).' - dataset.accelPosY.f_b__i_b_true(2,:).';
d((1:K) + 3*K) = dataset.accelNegY.f_b__i_b_meas(2,:).' - dataset.accelNegY.f_b__i_b_true(2,:).';

% Build from Z Tests
G((1:K) + 4*K, (1:4) + 2*4) = [ones(K, 1), dataset.accelPosZ.f_b__i_b_true.'];
G((1:K) + 5*K, (1:4) + 2*4) = [ones(K, 1), dataset.accelNegZ.f_b__i_b_true.'];
d((1:K) + 4*K) = dataset.accelPosZ.f_b__i_b_meas(3,:).' - dataset.accelPosZ.f_b__i_b_true(3,:).';
d((1:K) + 5*K) = dataset.accelNegZ.f_b__i_b_meas(3,:).' - dataset.accelNegZ.f_b__i_b_true(3,:).';


end

