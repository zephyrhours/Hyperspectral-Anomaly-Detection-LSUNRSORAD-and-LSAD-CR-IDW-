function result = Local_RX(Data, win_out, win_in)
%HYPERRX RX anomaly detector
%   hyperRxDetector performs the RX anomaly detector
%
% Usage
%   [result] = hyperRxDetector(Data, window, lambda)
% Inputs
%   Data - 3D data matrix (num_row x num_col x num_dim)
%   win_out - spatial size window of outer(e.g., 3, 5, 7, 9,...)
%   win_in - spatial size window of inner(e.g., 3, 5, 7, 9,...)
%   lambda - regularization parameter
% Outputs
%   result - Detector output (num_row x num_col)
%  

[a b c] = size(Data);
result = zeros(a, b);
t = fix(win_out/2);
t1 = fix(win_in/2);
M = win_out^2;

% padding avoid edges
DataTest = zeros(3*a, 3*b, c);
DataTest(a+1:2*a, b+1:2*b, :) = Data;
DataTest(a+1:2*a, 1:b, :) = Data(:, b:-1:1, :);
DataTest(a+1:2*a, 2*b+1:3*b, :) = Data(:, b:-1:1, :);
DataTest(1:a, :, :) = DataTest(2*a:-1:(a+1), :, :);
DataTest(2*a+1:3*a, :, :) = DataTest(2*a:-1:(a+1), :, :);

for i = 1+b: 2*b 
    for j = 1+a: 2*a
        block = DataTest(j-t: j+t, i-t: i+t, :);
        y = squeeze(DataTest(j, i, :)).';   % 1 x dim
        block(t-t1+1:t+t1+1, t-t1+1:t+t1+1, :) = NaN;
        block = reshape(block, M, c);
        block(isnan(block(:, 1)), :) = [];
        H = block';  % num_dim x num_sam
        Sigma = (H * H');
        Sigma_inv = pinv(Sigma); 
        result(j-a, i-b) = y * Sigma_inv * y';
    end
end