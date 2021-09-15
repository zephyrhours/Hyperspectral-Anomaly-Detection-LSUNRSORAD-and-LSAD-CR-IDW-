function RDet=fun_LSAD( hsi,win_size)
%% Compiled by Zephyr Hou on 2017-12-12
% 
% LSAD( Local Summation Anomaly Detection )
%
% The References:
% 1.¡¶A spectral-spatial based loacl summation anomaly detection 
%    method for hyspectral images¡·
%% Function Usage
% Input:
%      hsi£ºthe data cube, size of rows x columns x bands;
%      win_size £º the detection local windows size with the size
%              2*w+1(w=0,1,2,...);
% Output:
%      Res_Det:the result of LSADBCR Detection, size of  rows x columns
%% -----------------------LSAD Main Function ------------------------------ 
[rows,cols,bands]=size(hsi);

%%% nomalized
hsi1=reshape(hsi,rows*cols,bands);
max_y = max(max(hsi1));
min_y = min(min(hsi1));
hsi1 = (hsi1-min_y)/(max_y-min_y);
hsi=reshape(hsi1,rows,cols,bands);

%%% parameter
if nargin < 2
    win_size=3;
end

w=(win_size-1)/2;
% randomly choose the pixels to expand 2*w layers of pixels on the edge of image
Exp_hsi=zeros(2*w+rows,2*w+cols,bands);
for i=1:rows+4*w
    for j=1:cols+4*w
        ind_x=randperm(rows);
        ind_y=randperm(cols);
        Exp_hsi(i,j,:)=hsi(ind_x(1),ind_y(1),:);
    end
end
Exp_hsi(2*w+1:2*w+rows,2*w+1:2*w+cols,:)=hsi;
 
RDet=zeros(rows,cols);
R_Det1=zeros(win_size,win_size);

for i=2*w+1:2*w+rows
    for j=2*w+1:2*w+cols
        Img_block=Exp_hsi(i-2*w:i+2*w,j-2*w:j+2*w,:);
        % Img_block(2*w+1,2*w+1,:)=zeros(1,bands);
        for ki=1:win_size
            for kj=1:win_size
                Img_blocki=Img_block(ki:ki+win_size-1,kj:kj+win_size-1,:);
                re_Img_blocki=reshape(Img_blocki,win_size*win_size,bands);
                cov_blocki=cov(re_Img_blocki);
                % PCA
                [eig_XL,eig_Z]=eig(cov_blocki);
                [Deig_Z,ind]=sort(diag(eig_Z),'descend');
                D_eigXL=eig_XL(:,ind');
                inv_Cov=D_eigXL(:,1:5)*pinv(diag(Deig_Z(1:5)))*D_eigXL(:,1:5)';
                Y=reshape(Img_block(2*w+1,2*w+1,:),bands,1)-reshape(mean(re_Img_blocki),bands,1);           
                DetY=Y'*inv_Cov*Y;
                R_Det1(ki,kj)=DetY;
            end
        end
        RDet(i-2*w,j-2*w)=sum(sum(R_Det1));
    end
end
end

