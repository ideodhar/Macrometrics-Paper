clear all
close all
clc

nhor=30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BDM baseline output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load asym_results_BDM.mat
for hor_Msum=1:nhor


% Count # of points above 45degree line:
count=sum(squeeze(MsumS(:,2,hor_Msum))>squeeze(MsumS(:,1,hor_Msum)));

results_prob(hor_Msum)= 1- count/length(MsumS);
end

orig= results_prob;

%Modifying K=30

load asym_results_BDM_30lags.mat
for hor_Msum=1:nhor


% Count # of points above 45degree line:
count=sum(squeeze(MsumS(:,2,hor_Msum))>squeeze(MsumS(:,1,hor_Msum)));

results_prob(hor_Msum)= 1- count/length(MsumS);
end

orig_30= results_prob;

%Modifying K=20

load asym_results_BDM_20lags.mat
nhor=20;
for hor_Msum=1:nhor


% Count # of points above 45degree line:
count=sum(squeeze(MsumS(:,2,hor_Msum))>squeeze(MsumS(:,1,hor_Msum)));

results_prob(hor_Msum)= 1- count/length(MsumS);
end

orig_20= results_prob;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUTS from running over the full RZ sample, 1890-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load asym_results_fullRZsample.mat
for hor_Msum=1:nhor


% Count # of points above 45degree line:
count=sum(squeeze(MsumS(:,2,hor_Msum))>squeeze(MsumS(:,1,hor_Msum)));

results_prob(hor_Msum)= 1- count/length(MsumS);
end

origfullsample= results_prob;


load asym_results_fullRZsample_30lags.mat
for hor_Msum=1:nhor


% Count # of points above 45degree line:
count=sum(squeeze(MsumS(:,2,hor_Msum))>squeeze(MsumS(:,1,hor_Msum)));

results_prob(hor_Msum)= 1- count/length(MsumS);
end

orig_30fullsample= results_prob;


load asym_results_fullRZsample_20lags.mat
nhor=20;
for hor_Msum=1:nhor


% Count # of points above 45degree line:
count=sum(squeeze(MsumS(:,2,hor_Msum))>squeeze(MsumS(:,1,hor_Msum)));

results_prob(hor_Msum)= 1- count/length(MsumS);
end

orig_20fullsample= results_prob;



nhor=30;

figure(7)
plot(1:1:nhor, orig, 1:1:nhor, origfullsample, 'LineWidth', 1.5)
hold on
plot(1:1:nhor, 0.9*ones(nhor,1), 'k')
hold on 
plot(1:1:nhor, 0.95*ones(nhor,1), 'k')
xlabel('horizon, h')
title('P(msum neg>msum pos) at hor=h')
legend('BDM sample: 1901q2-2015q4', 'RZ full sample: 1890q1-2015q4', '', '')



figure(8)
subplot(1,2,1)
plot(1:1:nhor, orig, 'b', 1:1:nhor, orig_30 ,'r--', 'LineWidth', 2)
hold on
plot( orig_20, 'c-*', 'LineWidth', 2)
hold on
plot(1:1:nhor, 0.9*ones(nhor,1), 'k', 'LineWidth', 1)
hold on 
plot(1:1:nhor, 0.95*ones(nhor,1), 'k', 'LineWidth', 1)
xlabel('horizon, h')
title('1901q2-2015q4')

subplot(1,2,2)
plot(1:1:nhor, origfullsample, 'b', 1:1:nhor, orig_30fullsample , 'r--', 'LineWidth', 2)
hold on
plot( orig_20fullsample, 'c-*', 'LineWidth', 2)
hold on
plot(1:1:nhor, 0.9*ones(nhor,1), 'k', 'LineWidth', 1)
hold on 
plot(1:1:nhor, 0.95*ones(nhor,1), 'k', 'LineWidth', 1)
xlabel('horizon, h')
legend('Original specification w/ 45lags', 'w/ 30 lags', 'w/ 20 lags', '', '')
title('1890q1-2015q4')



