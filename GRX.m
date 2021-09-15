function    GRX_Detect= GRX( hsi )
% Complied by Zephyr Hou on 2017-03-01
%
% function usage
% Global RX
% GRX_Detect= GRX( hsi )
%
% Input:
%     hsi: the data cube, size of rows x columns x bands
% Output:
%     GRX_Detect: the result of Global RX Detection, size of rows x columns
% 

[rows,cols,bands]=size(hsi);
X=reshape(hsi,rows*cols,bands);

X_mean = mean(X)';   % the average value of each bands, bands x 1
cov_X=cov(X);       % the covariance value between bands, bands x bands
Y = X' - repmat(X_mean,[1,rows*cols]);
GRX_Detect=zeros(1,rows*cols);
cov_inv=pinv(cov_X);
for i=1:rows*cols
    GRX_Detect(1,i)=Y(:,i)'*cov_inv*Y(:,i);
end
GRX_Detect=reshape(GRX_Detect,rows,cols);
end

