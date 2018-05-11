% This script tries to run the three population model

clear all, close all, clc

vddata = load ('eight_weeks_data_2_18.m'); % load concatenated data set (long Mx2 matrix
% all weeks added vertically, first column is week number
naive_vddata = load('naive_data.m');
%% Remove Cohort 3
coh_ind = vddata(:, 4) == 3;
%%
vddata = vddata(coh_ind == 0, :);
% coh_ind8 = vddata(:,4) == 8;
% vddata = vddata(coh_ind8 == 0, :);
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

% %[new_vddata] = fit_by_var_by_time(vddata);
% %figure(10)
% hold off
% %plot(dose, new_vddata(:,5),'o');
% %ylim([0 100])
% xlabel('dose(uM)', 'FontSize', 20)
% ylabel('Weight by Variance','FontSize', 20)
% title('Plot of Proposed Weighting Scheme', 'FontSize', 20)

%% Three Population Model fractions shift in time

options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
%define bounds for all parameters
            % m1 cen1 m2 cen 2 wk1sens wk1res etc...
paramslb = zeros([22 1]);
paramsub = [ Inf; Inf; Inf; Inf; Inf; Inf; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1];

% sets up initial values of beta2new
params0 = [ .1; 17; .04; 35; .04; 65; .3; .3; .3; .3; .3; .3; .3; .3; .3; .3; .3; .3; .3; .3; .3; .3];
%  res2sq=ones([n 1]);
[beta3new, resnorm3, residuals3] = lsqnonlin(@fit_simp3pop8weeks,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    dose,...
    var,...
    wknum,...
    Vmaxall);
coeffs3 = beta3new(1:6);


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
v_model3 = model3popallweeksnormed( dose, beta3new, Vmaxbyweek, nsize);
sigma3 = std(residuals3)
error3 = sum(abs(residuals3))./n
RSS3 = sum(residuals3.^2);
p3 = 22;
AIC3 = n*log(RSS3./n) + 2*p3
BIC3 = -n*log(RSS3./n) + p3*log(n)
chi_squared3 = sum((residuals3.^2)./v_model3)


DOF3 = n-p3;
MSE3 = chi_squared3./ DOF3
F_stat1v3alt = ((RSS1-RSS3)./(DOF1-DOF3))./(RSS3./(DOF3))
p1v3alt = 1-fcdf(F_stat1v3alt, DOF1-DOF3, DOF3)

% Comparison of one population changing in time with two population model
F_statLD50v3alt = ((RSSLD50-RSS3)./(DOFLD50-DOF3))./(RSS3./(DOF3))
p3vLD50alt = 1-fcdf(F_statLD50v3alt, DOFLD50-DOF3, DOF3)


F_stat2v3alt = ((RSS2-RSS3)./(DOF2-DOF3))./(RSS3./(DOF3))
p2v3alt = 1-fcdf(F_stat2v3alt, DOF2-DOF3, DOF3)


hold off
plot(dose, var, 'bo')
hold on
plot(dose, v_model3, 'ro')

figure(11)
hold off
plot3(dose, vddata(:,1), v_model3, 'ro', 'LineWidth',5)
hold on
plot3(dose, vddata(:,1), var, 'bo')
xlabel(' Dose Doxorubicin (uM)', 'FontSize', 12)
zlabel('Model Prediction and Experiment Viability', 'FontSize',12)
ylabel('Time (Weeks Post Treatment)','FontSize', 12)
%%
for i = 1:8
    fsens(1) = fsens0;
fsens(i+1) = beta3new(6+(2.*i-1));
end

for i = 1:8
    ftol(1)=ftol0;
    ftol(i+1)= beta3new(6+(2.*i));
end

for i = 1:8
    fres(1) = fres0;
    fres(i+1) = 1- ftol(i+1)-fsens(i+1);
end
time = 0:1:8;
%%
beta3new8 = beta3new;
[lowerlim3boot, upperlim3boot] = BSerrorinparams3(residuals3, beta3new, vddata);
%%
[lowerlim3bootnaive, upperlim3bootnaive] = BSerrorinparams3naive(residuals3naive, dosenc, beta3naive, naive_vddata, beta3new, Vmaxnaiveavg);
%%
errorbarlength3boot = upperlim3boot - lowerlim3boot

%%
figure(11)
hold off
errorbar(time(2:9), fsens(2:9), errorbarlength3boot(7:2:21)./2, 'yo-', 'LineWidth', 1.5)
hold on
errorbar(time(2:9), ftol(2:9), errorbarlength3boot(8:2:22)./2, 'go-', 'LineWidth', 1.5)
errorbar(time(2:9), fres(2:9), (errorbarlength3boot(8:2:22) + errorbarlength3boot(7:2:21))./4, 'bo-', 'LineWidth', 1.5)
%plot(time, fsens, '-yo', 'LineWidth', 2.5)
hold on
%plot(time, ftol, '-go', 'LineWidth', 2.5)
%plot(time, fres, '-bo', 'LineWidth', 2.5)
xlabel ('Time (Weeks Post Treatment)', 'FontSize', 16)
ylabel ('Fraction of Cells', 'FontSize', 16)
title (' Three Population Model over Eight Weeks', 'FontSize', 16)
legend ('Sensitive Fraction LD50 = 19.6 uM', 'Tolerant Fraction LD50 = 69.1 uM',' Resistant Fraction LD50 = 133.6 uM')
set(gca,'LineWidth',1.5,'FontSize',16)