%% Remote Sensing£¬Compiled by ZephyrHou on 2018-11-24
clc;clear;close all;
%% load dataset

%%% open the image dataset
load('Salians_syn.mat')
load('Salians_syn_gt.mat')

[rows,cols,bands]=size(hsi);
%% compared other methods
tic
R1=GRX(hsi);
t1=toc;
disp(['GRX Running Time£º',num2str(t1)])

tic
R2=Local_RX(hsi,5,3);
t2=toc;
disp(['LRX Running Time£º',num2str(t2)])

tic
R3=UNRS_ID_Detect(hsi,5,3,0.01);
t3=toc;
disp(['UNRS Running Time£º',num2str(t3)])

tic
R4=CRD_Detect(hsi,5,3,0.01);
t4=toc;
disp(['CRD Running Time£º',num2str(t4)])

tic
R5=fun_LSAD(hsi,5);
t5=toc;
disp(['LSAD Running Time£º',num2str(t5)])


%% Proposed method

tic
R6=fun_LSUNRS_ED_ORAD(hsi,5,3,0.01);
t6=toc;
disp(['LSUNRSORAD Running Time£º',num2str(t6)])

tic
R7=fun_LSAD_CR_IDW(hsi,5,3,0.01);
t7=toc;
disp(['LSAD-CR-IDW Running Time£º',num2str(t7)])

%% Show the detected results
label_value=reshape(hsi_gt,1,rows*cols);

R1value = reshape(R1,1,rows*cols);
R2value = reshape(R2,1,rows*cols);
R3value = reshape(R3,1,rows*cols);
R4value = reshape(R4,1,rows*cols);
R5value = reshape(R5,1,rows*cols);
R6value = reshape(R6,1,rows*cols);
R7value = reshape(R7,1,rows*cols);

[FA1,PD1] = perfcurve(label_value,R1value,'1') ;
[FA2,PD2] = perfcurve(label_value,R2value,'1') ;
[FA3,PD3] = perfcurve(label_value,R3value,'1') ;
[FA4,PD4] = perfcurve(label_value,R4value,'1') ;
[FA5,PD5] = perfcurve(label_value,R5value,'1') ;
[FA6,PD6] = perfcurve(label_value,R6value,'1') ;
[FA7,PD7] = perfcurve(label_value,R7value,'1') ;


%-------------------------------------------------------------------------

AUC1=-sum((FA1(1:end-1)-FA1(2:end)).*(PD1(2:end)+PD1(1:end-1))/2);
AUC2=-sum((FA2(1:end-1)-FA2(2:end)).*(PD2(2:end)+PD2(1:end-1))/2);
AUC3=-sum((FA3(1:end-1)-FA3(2:end)).*(PD3(2:end)+PD3(1:end-1))/2);
AUC4=-sum((FA4(1:end-1)-FA4(2:end)).*(PD4(2:end)+PD4(1:end-1))/2);
AUC5=-sum((FA5(1:end-1)-FA5(2:end)).*(PD5(2:end)+PD5(1:end-1))/2);
AUC6=-sum((FA6(1:end-1)-FA6(2:end)).*(PD6(2:end)+PD6(1:end-1))/2);
AUC7=-sum((FA7(1:end-1)-FA7(2:end)).*(PD7(2:end)+PD7(1:end-1))/2);

disp('-------------------------------------------------------------------')
disp('AUC Values and Running Times:')
disp('GRX')
disp(['AUC:     ',num2str(AUC1),'          Time:     ',num2str(t1)])
disp('LRX')
disp(['AUC:     ',num2str(AUC2),'          Time:     ',num2str(t2)])
disp('UNRS')
disp(['AUC:     ',num2str(AUC3),'          Time:     ',num2str(t3)])
disp('CRD')
disp(['AUC:     ',num2str(AUC4),'          Time:     ',num2str(t4)])
disp('LSAD')
disp(['AUC:     ',num2str(AUC5),'          Time:     ',num2str(t5)])
disp('LSUNRSORAD')
disp(['AUC:     ',num2str(AUC6),'          Time:     ',num2str(t6)])
disp('LSAD-CR-IDW')
disp(['AUC:     ',num2str(AUC7),'          Time:     ',num2str(t7)])
disp('-------------------------------------------------------------------')

%%
figure(1);
semilogx(FA1, PD1, 'c-', 'LineWidth', 2);  hold on
semilogx(FA2, PD2, 'm-', 'LineWidth', 2);  hold on
semilogx(FA3, PD3, 'g-', 'LineWidth', 2);  hold on
semilogx(FA4, PD4, 'y-', 'LineWidth', 2);  hold on
semilogx(FA5, PD5, 'b-', 'LineWidth', 2);  hold on
semilogx(FA6, PD6, 'r-', 'LineWidth', 2);  hold on
semilogx(FA7, PD7, 'k-', 'LineWidth', 2);  hold on

axis([0.65e-4,1,0,1])
xlabel('False alarm rate'); ylabel('Probability of detection');
legend('GRX','LRX','UNRS','CRD','LSAD','LSUNRSORAD','LSAD-CR-IDW','location','southeast')
title('Salinas');

%%
figure(2);
plot(FA1, PD1, 'c-', 'LineWidth', 2);  hold on
plot(FA2, PD2, 'm-', 'LineWidth', 2);  hold on
plot(FA3, PD3, 'g-', 'LineWidth', 2);  hold on
plot(FA4, PD4, 'y-', 'LineWidth', 2);  hold on
plot(FA5, PD5, 'b-', 'LineWidth', 2);  hold on
plot(FA6, PD6, 'r-', 'LineWidth', 2);  hold on
plot(FA7, PD7, 'k-', 'LineWidth', 2);  hold on

axis([0.65e-4,1,0,1])
xlabel('False alarm rate'); ylabel('Probability of detection');
legend('GRX','LRX','UNRS','CRD','LSAD','LSUNRSORAD','LSAD-CR-IDW','location','southeast')
title('Salinas');






