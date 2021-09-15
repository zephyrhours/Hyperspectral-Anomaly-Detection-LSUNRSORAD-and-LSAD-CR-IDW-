function RDet=fun_LSAD_CR_IDW(hsi,win_out,win_in,lambda )
%% Compiled by Zephyr Hou on 2017-12-13
%
% DWLSADBCR_DW(Dual Windows Local Summation Anomaly Detection Based on Collaborative  
% Representation with Distance-weighted)
%
% The References:
% 1.¡¶Hyperspectral Anomaly Detection Using Outlier Removal
%    from Collaborative Representation¡·
% 2.¡¶A spectral-spatial based loacl summation anomaly detection 
%    method for hyspectral images¡·
%% Function Usage
% Input:
%      hsi£ºthe data cube, size of rows x columns x bands;
%      win_out £º the detection local outer windows size with the size
%              2*w+1(w=0,1,2,...);
%      win_in £º the detection local inner windows size with the size
%              2*w1+1(w=0,1,2,...);
%      lambda£ºa parameter to adjust the norm of weight vectors in the
%               image depending on its backgrounds;
% Output:
%      RDet:the result of DWLASDBCR-DW Detection, size of  rows x columns
%% ---------------------DWLSADBCR-DW ³ÌÐò----------------------------
[rows,cols,bands]=size(hsi);
%%% Normalization
hsi1=reshape(hsi,rows*cols,bands);
maxVal = max(max(hsi1));
minVal = min(min(hsi1));
hsi1 = (hsi1-minVal)/(maxVal-minVal);
hsi=reshape(hsi1,rows,cols,bands);

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
Taoy=zeros(num_S,num_S);
%% ================================================
% IDW 
M = win_out^2;
tempWinX=zeros(win_out,win_out);
tempWinY=zeros(win_out,win_out);
for i=1:win_out
    tempWinX(i,:)=-t:t;
    tempWinY(:,i)=t:-1:-t;
end

Dd=sqrt(tempWinX.^2+tempWinY.^2);
Dd(t-t1+1:t+t1+1, t-t1+1:t+t1+1) = NaN;
IDW = reshape(Dd, M, 1);
IDW(isnan(IDW(:,1)),:) = [];
SumW=sum(IDW.^-2);
IDW=IDW.^-2/SumW;

%% ===============================================================

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
                H = block;  % num_sam x num_dim        
                temp=CenPix-H';
                for ii=1:num_S
                    Taoy(ii,ii)=norm(temp(:,ii)); % Tikhonov
                end
                Taoy=IDW.*Taoy;
%                 Taoy=SDW.*Taoy;
                W=pinv(H*H'+lambda*(Taoy'*Taoy))*H*CenPix;  

                DetY= norm(CenPix-H'*W);
                tempR(ki+t1+1,kj+t1+1)=DetY;

            end
        end
    
        RDet(i-rows,j-cols)=sum(sum(tempR)) ;                 
    end
end
t=toc;
disp(['LSAD-CR-IDW running time£º',num2str(t)]);
end


