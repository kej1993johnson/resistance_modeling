% % Script plots histograms of the simulated parameter sets
% paramset2 = load('beta2newsim.mat');
% paramset2 = struct2cell(paramset2);
% paramset2 = cell2mat(paramset2);
paramset2 = beta2newsim;

% paramset3 = load('beta3newsim.mat');
% paramset3 = struct2cell(paramset3);
% paramset3 = cell2mat(paramset3);
% beta3newsimnorm = paramset3;
% 
% paramset3norm = load('beta3newnormsim.mat');
% paramset3norm = struct2cell(paramset3norm);
% paramset3norm = cell2mat(paramset3norm);
%%


figure(1)
hold off
subplot (1,2,1) % plots beta2new sim centers
hist(paramset2(2,:),100);
xlabel('Possible sensitive LD50 parameter fits','FontSize', 18)
ylabel('Frequency','FontSize', 18)



subplot (1,2,2) 
hold off
hist(paramset2(4,:),100);
xlabel('Possible resistant LD50 parameter fits', 'FontSize', 18)
ylabel ('Frequency', 'FontSize', 18)

%%


hold off
figure(2)
subplot(2,6,1)
hist(paramset2(5,:),100);
xlabel('Week 1 Sensitive','FontSize', 14)
xlim([0 1])
ylabel('Frequency')
% legend('mean wk1 fraction sens = 0.79')

subplot(2,6,2)
hist(paramset2(6,:),100);
xlabel('Week 2 Sensitive','FontSize', 14)
xlim([0 1])
ylabel('Frequency')
%legend('mean wk2 fraction sens = 0.28')

subplot(2,6,3)
hist(paramset2(7,:),100);
xlabel('Week 3 Sensitive','FontSize', 14)
xlim([0 1])
ylabel('Frequency')
% legend('mean wk3 fraction sens = 0.20')

subplot(2,6,4)
hist(paramset2(8,:),100);
xlabel('Week 4 Sensitive','FontSize', 14)
xlim([0 1])
ylabel('Frequency')
%legend('mean wk4 fraction sens = 0.43')

subplot(2,6,5)
hist(paramset2(9,:),100);
xlabel('Week 5 Sensitive','FontSize', 14)
xlim([0 1])
ylabel('Frequency')
%legend('mean wk5 fraction sens = 0.62')

subplot(2,6,6)
hist(paramset2(10,:),100);
xlabel('Week 6 Sensitive','FontSize', 14)
xlim([0 1])
ylabel('Frequency')
%legend('mean wk6 fraction sens = 0.66')

subplot(2,6,7)
hist(paramset2(11,:),100);
xlabel('Week 7 Sensitive','FontSize', 14)
xlim([0 1])
ylabel('Frequency')
%legend('mean wk7 fraction sens = 0.67')

subplot(2,6,8)
hist(paramset2(12,:),100);
xlabel('Week 8 Sensitive','FontSize', 14)
xlim([0 1])
ylabel('Frequency')
%legend('mean wk8 fraction sens = 0.62')

subplot(2,6,9)
hist(paramset2(13,:),100);
xlabel('Week 9 Sensitive','FontSize', 14)
xlim([0 1])
ylabel('Frequency')
%legend('mean wk9 fraction sens = 0.61')

subplot(2,6,10)
hist(paramset2(14,:),100);
xlabel('Week 10 Sensitive','FontSize', 14)
xlim([0 1])
ylabel('Frequency')
%legend('mean wk10 fraction sens = 0.69')

subplot(2,6,11)
hist(paramset2(15,:),100);
xlabel('Week 11 Sensitive','FontSize', 14)
xlim([0 1])
ylabel('Frequency')
%legend('mean wk11 fraction sens = 0.57')

subplot(2,6,12)
hist(paramset2(16,:),100);
xlabel('Week 12 Sensitive','FontSize', 14)
xlim([0 1])
ylabel('Frequency')
%legend('mean wk12 fraction sens = 0.58')
%%
hold off
figure (3)
subplot(1,3,1)
hist(paramset2(7,:),100);
xlabel('Week 3 Sensitive','FontSize', 14)
xlim([0 1])
ylabel('Frequency')
legend('mean = 0.42')

subplot(1,3,2)
hist(paramset2(10,:),100);
xlabel('Week 6 Sensitive','FontSize', 14)
xlim([0 1])
ylabel('Frequency')
legend('mean = 0.52')

subplot(1,3,3)
hist(paramset2(14,:),100);
xlabel('Week 10 Sensitive','FontSize', 14)
xlim([0 1])
ylabel('Frequency')
legend('mean = 0.64')


%% Do the same for the three population model


figure(4)
hold off
subplot (1,3,1) % plots beta3new sim centers
hist(paramset3(2,:),100);
xlabel('Possible sensitive LD50 parameter fits')
ylabel('Frequency')
%legend('mean LD50 sens = 31.5 uM')


subplot (1,3,2) 
hist(paramset3(4,:),100);
xlabel('Possible tolerant LD50 parameter fits')
ylabel ('Frequency')
%legend('mean LD50 res = 65.0 uM')

subplot (1,3,3) 
hist(paramset3(6,:),100);
xlabel('Possible resistant LD50 parameter fits')
ylabel ('Frequency')
%legend('mean LD50 res = 65.0 uM')

%%

hold off
figure(5)
subplot(2,6,1)
hist(paramset3(8,:),100);
xlabel('Week 1')
ylabel('Frequency')
%legend('mean wk1 fraction tol = 0.79')

subplot(2,6,2)
hist(paramset3(10,:),100);
xlabel('Week 2')
ylabel('Frequency')
%legend('mean wk2 fraction tol = 0.28')

subplot(2,6,3)
hist(paramset3(12,:),100);
xlabel('Week 3')
ylabel('Frequency')
%legend('mean wk3 fraction tol = 0.20')

subplot(2,6,4)
hist(paramset3(14,:),100);
xlabel('Week 4')
ylabel('Frequency')
%legend('mean wk4 fraction tol = 0.43')

subplot(2,6,5)
hist(paramset3(16,:),100);
xlabel('Week 5')
ylabel('Frequency')
%legend('mean wk5 fraction tol = 0.62')

subplot(2,6,6)
hist(paramset3(18,:),100);
xlabel('Week 6')
ylabel('Frequency')
%legend('mean wk6 fraction tol = 0.66')

subplot(2,6,7)
hist(paramset3(20,:),100);
xlabel('Week 7')
ylabel('Frequency')
%legend('mean wk7 fraction tol = 0.67')

subplot(2,6,8)
hist(paramset3(22,:),100);
xlabel('Week 8')
ylabel('Frequency')
%legend('mean wk8 fraction tol = 0.62')

subplot(2,6,9)
hist(paramset3(24,:),100);
xlabel('Week 9')
ylabel('Frequency')
%legend('mean wk9 fraction tol = 0.61')

subplot(2,6,10)
hist(paramset3(26,:),100);
xlabel('Week 10')
ylabel('Frequency')
%legend('mean wk10 fraction tol = 0.69')

subplot(2,6,11)
hist(paramset3(28,:),100);
xlabel('Week 11')
ylabel('Frequency')
%legend('mean wk11 fraction sens = 0.57')

subplot(2,6,12)
hist(paramset3(30,:),100);
xlabel('Week 12')
ylabel('Frequency')
%legend('mean wk12 tolerant sens = 0.58')

%%
hold off
figure (6)
subplot(1,3,1)
hist(paramset3(12,:),100);
xlabel('Week 3 Tolerant','FontSize', 14)
xlim([0 1])
ylabel('Frequency')
legend('mean = 0.08')

subplot(1,3,2)
hist(paramset3(18,:),100);
xlabel('Week 6 Tolerant','FontSize', 14)
xlim([0 1])
ylabel('Frequency')
legend('mean = 0.21')

subplot(1,3,3)
hist(paramset3(26,:),100);
xlabel('Week 10 Tolerant','FontSize', 14)
xlim([0 1])
ylabel('Frequency')
legend('mean = 0.32')

