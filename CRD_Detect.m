function result = CRD_Detect(Data, win_out, win_in, lambda)
%
% Compiled by Zephyr Hou on 2017-02-26
% Collaborative Representation for Hyperspectral Anomaly Detector
%
% Refering to the paper of Collaborative Representation for Hyperspectral Anomaly Detector
% 
% Usage
%   [result] = CRD_Detect(Data, window, lambda)
% Inputs
%   Data - 3D data matrix (num_row x num_col x num_dim)
%   window - spatial size window (e.g., 3, 5, 7, 9,...)
%   lambda - regularization parameter
% Outputs
%   result - Detector output (num_row x num_col)
%  

[a b c] = size(Data);        
result = zeros(a, b);
t = fix(win_out/2);
t1 = fix(win_in/2);
M = win_out^2;
num_sam=win_out*win_out-win_in*win_in;

% % padding avoid edges
% DataTest = zeros(3*a, 3*b, c);
% DataTest(a+1:2*a, b+1:2*b, :) = Data;
% DataTest(a+1:2*a, 1:b, :) = Data(:, b:-1:1, :);
% DataTest(a+1:2*a, 2*b+1:3*b, :) = Data(:, b:-1:1, :);
% DataTest(1:a, :, :) = DataTest(2*a:-1:(a+1), :, :);
% DataTest(2*a+1:3*a, :, :) = DataTest(2*a:-1:(a+1), :, :);

% padding zeros to avoid edges
DataTest = zeros(3*a, 3*b, c);
DataTest(a+1:2*a, b+1:2*b, :) = Data;
DataTest(a+1:2*a, 1:b, :) = zeros(a,b,c);
DataTest(a+1:2*a, 2*b+1:3*b, :) = zeros(a,b,c);
DataTest(1:a,:,:)=zeros(a,3*b,c);
DataTest(2*a+1:3*a,:,:)=zeros(a,3*b,c);



for i = 1+b: 2*b 
    for j = 1+a: 2*a
        block = DataTest(j-t: j+t, i-t: i+t, :);
        y = squeeze(DataTest(j, i, :)).';% 1 x num_dim
        block(t-t1+1:t+t1+1, t-t1+1:t+t1+1, :) = NaN;
        block = reshape(block, M, c);
        block(isnan(block(:, 1)), :) = [];
        Xs = block';  % num_dim x num_sam
        
        weights=pinv(Xs'*Xs+lambda*eye(num_sam))*Xs'*y';% num_sam x 1(formula 3)¡¾CRD¡¿ 
        
        y_hat = (Xs*weights(:))';  % 1 x num_dim
        result(j-a, i-b) = norm(y - y_hat, 2);
    end
end

