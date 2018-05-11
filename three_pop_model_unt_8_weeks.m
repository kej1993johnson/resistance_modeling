% This script tries to run the three population model

clear all, close all, clc

vddata = load ('eight_weeks_data_2_18_untreated.m'); % load concatenated data set (long Mx2 matrix
% all weeks added vertically, first column is week number
naive_vddata = load('naive_data.m');

%% Allocate data

dose = vddata (:,2);
var = vddata (:,3);
cohort_number = vddata(:,4);
n = length(dose); % finds number of data points total
wk = vddata(:,1);
[ Vmaxbyweek, Vmaxweekavg, ninweek, wknum, Vmaxall] = findVmaxandsize(vddata);
% Function tells you the average Max Viability per week (Vmaxbyweek), the
% average over all of the data (Vmaxweekavg), the number of dose-response
% curves in each week (so number of different cohorts in that week) the
% number of actual data points in each week of data (wknum(:,2), generally
% 12 x ninweek but some weeks have less data) and Vmaxall which puts the
% Vmax for each week (Vmax by week) alongside each dose and viability point
% to be used to calculate the fsens parameter for that week


% do the same for the naive data
dosenc = naive_vddata(:,2);
varnc = naive_vddata(:,3);

dose0ind = naive_vddata(:,2) == 0;
 Vmaxcohortnaivedata = naive_vddata(dose0ind,:);
 Vmaxnaiveavg = mean(Vmaxcohortnaivedata(:,3));
 
nsize = wknum(:,2); % 12 x 1 matrix of number of data points in each week
cohorts = unique(vddata(:,4));% gives cohort numbers (3,4,6,8,9)
num_cohorts = length(cohorts); % tells number of cohorts
[kstd, Wstd, std_dev_model_fun] = fittechnicalerrors( dose);

[new_vddata] = fit_by_var_by_time(vddata);
figure(10)
hold off
plot(dose, new_vddata(:,5),'o');
%ylim([0 100])
xlabel('dose(uM)', 'FontSize', 20)
ylabel('Weight by Variance','FontSize', 20)
title('Plot of Proposed Weighting Scheme', 'FontSize', 20)

%% Three Population Model fractions shift in time

options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
%define bounds for all parameters
            % m1 cen1 m2 cen 2 wk1sens wk1res etc...
paramslb = zeros([16 1]);
paramsub = [ 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1;];

% sets up initial values of beta2new
params0 = [.3; .3; .3; .3; .3; .3; .3; .3; .3; .3; .3; .3; .3; .3; .3; .3];
%  res2sq=ones([n 1]);
coeffs3 = load('beta3new8.mat');
coeffs3 = struct2cell(coeffs3);
coeffs3 = cell2mat(coeffs3);
[beta3unt, resnorm3, residuals3] = lsqnonlin(@fit_simp3pop8weeks_untreated,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    dose,...
    var,...
    wknum,...
    Vmaxall,...
    coeffs3);

%% Find Naive three pop
params0 = [.7,.1];
paramslb = [0,0];
paramsub = [1,1];

[beta3naive, resnorm3naive, residuals3naive] = lsqnonlin(@fit_simp3pop6weeksnaive,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    dosenc,...
    varnc,...
    coeffs3);

fsens0 = beta3naive(1);
fres0 = beta3naive(2);
ftol0 = 1-fsens0-fres0;
%% Compute and Output Two Population Statistics and Model Comparison

sigma3 = std(residuals3)
error3 = sum(abs(residuals3))./n
RSS3 = sum(residuals3.^2);
p3 = 22;
AIC3 = n*log(RSS3./n) + 2*p3
BIC3 = -n*log(RSS3./n) + p3*log(n)

DOF3 = n-p3;

for i = 1:8
    fsensunt(1) = fsens0;
fsensunt(i+1) = beta3unt((2.*i-1));
end

for i = 1:8
    fresunt(1)=fres0;
    fresunt(i+1)= beta3unt((2.*i));
end

for i = 1:8
    ftolunt(1) = ftol0;
    ftolunt(i+1) = 1- fresunt(i+1)-fsensunt(i+1);
end
time = 0:1:8;
%%

figure(11)
hold off
plot(time, fsensunt, 'yo-', 'LineWidth', 1.5)
hold on
plot(time, ftolunt, 'go-', 'LineWidth', 1.5)
plot(time, fresunt, 'bo-', 'LineWidth', 1.5)
ylim([-.2 1])
xlabel ('Time (Weeks Post Treatment)', 'FontSize', 16)
ylabel ('Fraction of Cells', 'FontSize', 16)
title (' Three Population Model Untreated Cells First Six Weeks', 'FontSize', 16)
legend ('Sensitive Fraction LD50 = 20.3 uM', 'Tolerant Fraction LD50 = 70.9 uM',' Resistant Fraction LD50 = 133.0 uM')
set(gca,'LineWidth',1.5,'FontSize',16)