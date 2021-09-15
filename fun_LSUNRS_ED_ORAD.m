function RDet = fun_LSUNRS_ED_ORAD(hsi,win_out,win_in,lambda )
%% Combiled by Zephyr Hou on 2018-11-10
% Reference£º
% ¡¶Hyperspectral Anomaly Detection Using Outlier Removal from Collaborative Representation¡·
% ¡¶Unsupervised Nearest Regularized Subspace Anomaly Detection in Hyspectral Imagery¡·
% ¡¶A spectral-spatial based local summation anomaly detection method for hyspectral images¡·
%% Function Usage
% Inputs
%   hsi - 3D data matrix (num_row x num_col x num_dim)
%   win_out - spatial size of outer window (e.g., 3, 5, 7, 9,...)
%   win_in- spatial size of inner window (e.g., 3, 5, 7, 9,...)
%   lambda - regularization parameter
% Outputs
%   Rdet - Detect Map output (num_row x num_col)
%%  
[rows,cols,bands] = size(hsi);
%%% Normalization
hsi1=reshape(hsi,rows*cols,bands);
max_y = max(max(hsi1));
min_y = min(min(hsi1));
hsi1 = (hsi1-min_y)/(max_y-min_y);
hsi=reshape(hsi1,rows,cols,bands);
%%
% % parameter
% win_out=7;
% win_in=5;
% lambda=0.1;

%%% boundary expansion
t = fix(win_out/2);
t1 = fix(win_in/2);
M=win_out*win_out;

DataTest = zeros(3*rows, 3*cols, bands);
DataTest(rows+1:2*rows, cols+1:2*cols, :) = hsi;
DataTest(rows+1:2*rows, 1:cols, :) = hsi(:,cols:-1:1, :);
DataTest(rows+1:2*rows, 2*cols+1:3*cols, :) = hsi(:, cols:-1:1, :);
DataTest(1:rows, :, :) = DataTest(2*rows:-1:(rows+1), :, :);
DataTest(2*rows+1:3*rows, :, :) = DataTest(2*rows:-1:(rows+1), :, :);


%%
RDet=zeros(rows,cols);
tempR=zeros(win_in,win_in);
num_S=win_out*win_out-win_in*win_in;    

tic;
for i = 1+rows: 2*rows      
    for j = 1+cols: 2*cols  
        CenPix= squeeze(DataTest(i, j, :));   % dim x 1    
        for ki=-t1:t1
            for kj=-t1:t1
                block = DataTest(i+ki-t: i+ki+t, j+kj-t: j+kj+t, :);
                block(t-t1+1:t+t1+1, t-t1+1:t+t1+1, :) = NaN;
                block = reshape(block, M, bands);
                block(isnan(block(:, 1)), :) = [];   
                H = block';  % num_dim x num_sam
                temp1=sum(H);
                miu=mean(temp1);
                sigma=std(temp1);
                threshold_max=miu+2*sigma;
                threshold_min=miu-2*sigma;
                for ii=1:size(H,2)
                    if temp1(ii) < threshold_min || temp1(ii) > threshold_max
                        temp1(ii)=NaN;
                    end
                end
                H(:,isnan(temp1(:))) = [];
                num_sam=size(H,2);
                
                z = H - repmat(CenPix,1,size(H,2)); % num_dim x num_sam
                Gramm=z'*z;     % num_sam x num_sam
                
                tau=[];
                for k=1:num_sam
                    tau(1,k)=norm(z(:,k),2);     
                end
                
                tau_y=diag(tau); % num_sam x num_sam    (formula 14) 

                C=Gramm+lambda*(tau_y'*tau_y);
                
                invC=pinv(C);   % num_sam x num_sam
                W = sum(invC,2)/sum(sum(invC));     % num_sam *1
                DetY= norm(CenPix-H*W);
                tempR(ki+t1+1,kj+t1+1)=DetY;  
      
            end
        end
        RDet(i-rows,j-cols)=sum(sum(tempR)) ;                 
    end
end
t=toc;
% disp(['LSUNRS-ED-ORAD running time£º',num2str(t)]);

end

