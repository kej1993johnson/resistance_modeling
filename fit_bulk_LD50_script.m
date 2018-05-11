% This script takes data for all weeks and computes the average LD50s for
% each time point.  It plots the LD50s as a function of time. Then it does
% the same for each individual cohort, and plots the LD50s as a function of
% time on one graph

clear all, close all, clc


vddata = load ('allweeksdata11_22.m'); % load concatenated data set (long Mx2 matrix
% all weeks added vertically, first column is week number
naive_vddata = load('naive_data.m');

%%

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

nsize = wknum(:,2); % 12 x 1 matrix of number of data points in each week
cohorts = unique(vddata(:,4));% gives cohort numbers (3,4,6,8,9)
num_cohorts = length(cohorts); % tells number of cohorts
%% 12 weeks of data

options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
%define bounds for all parameters
            % m1 cen1 m2 cen 2 wk1sens wk1res etc...
paramslb = zeros([24 1]);
paramsub = Inf([24 1]);


params0 = [ .1; .1; .1; .1; .1; .1; .1; .1; .1; .1; .1; .1; 30; 30; 30; 30; 30; 30; 30; 30; 30; 30; 30; 30];
% outputs betaLD50 with first 12 paramters slope values, second 12 are
% center values in week order
[betaLD50, resnormLD50, residualsLD50s] = lsqnonlin(@fit_simpLD50,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    dose,...
    var,...
    wknum,...
    Vmaxall);

%% Find LD50 overall for naive data

dosen = naive_vddata (:,2);
varn = naive_vddata (:,3);
cohort_numbern = naive_vddata(:,4);
n_n = length(dosen); % finds number of data points total
wkn = naive_vddata(:,1);


options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
%define bounds for all parameters
            % m1 cen1 m2 cen 2 wk1sens wk1res etc...
paramslbn = zeros([2 1]);
paramsubn = Inf([2 1]);

dose0ind = naive_vddata(:,2) == 0;
 Vmaxnaivedata = naive_vddata(dose0ind,:);
 Vmaxnaiveavg = mean(Vmaxnaivedata(:,3));
params0n = [ .1; 30];
% outputs betaLD50 with first 12 paramters slope values, second 12 are
% center values in week order
[betaLD50naiveall, resnormLD50naive, residualsLD50snaive] = lsqnonlin(@fit_simp1popnaive,...
    params0n,...
    paramslbn,...
    paramsubn,...
    options,...
    dosen,...
    varn,...
    Vmaxnaiveavg);



%% Plot of LD50s

time(1, :) = 1:12;
figure (1)
% Plot of LD50s over time
plot(time, betaLD50(13:24), 'o', 'LineWidth', 4)
hold on
plot (0, betaLD50naiveall(2), 'go','LineWidth',4)
xlabel ('Time Post Treatment', 'FontSize', 20)
ylabel ('Bulk Population LD50', 'FontSize', 20)
ylim ([20,80])
title('Bulk Population Dynamics for All Data','FontSize', 20)

%% Model Fit

bulk_model_LD50 = model_bulk_LD50( dose, betaLD50, Vmaxall, wk );

figure(2)
% Plot of bulk model fits versus experiment
plot (dose, bulk_model_LD50,'o')
xlabel(' Dose Doxorubicin (uM)', 'FontSize', 20)
ylabel('Model Prediction of Viability', 'FontSize',20)

%% Compute All Statistics

actualresidualsLD50 = bulk_model_LD50-var;
sigmaLD50 = std(actualresidualsLD50)
errorLD50 = sum(abs(actualresidualsLD50))./n
chi_squaredLD50 = sum((actualresidualsLD50.^2)./bulk_model_LD50)
RSSLD50 = sum(actualresidualsLD50.^2);
pLD50 = 24;
AICLD50 = n*log(RSSLD50./n) + 2*pLD50
BICLD50 = -n*log(RSSLD50./n) + pLD50*log(n)
DOFLD50 = n-pLD50;
MSELD50 = chi_squaredLD50./ DOFLD50

%% Now find LD50s of the individual cohorts

for k=1:num_cohorts
    data_maniptest = vddata;
    index1 = data_maniptest(:,4) == cohorts(k);
    testingset = data_maniptest(index1,:);
    cell_of_testingsets{k} = testingset;   %makes separate cell of each cohort         
end   

for k = 1:num_cohorts
    testing_vddata = cell_of_testingsets(k);
    testing_vddata = cell2mat(testing_vddata);
    
 [ VmaxbyweekCVt, VmaxweekavgCVt, ninweekCVt, wknumCVt, VmaxallCVt] = findVmaxandsizetest(testing_vddata);
nsizeCVt = wknumCVt(:,2);

VmaxbyweekCVtcell{k} = VmaxbyweekCVt;
VmaxweekavgCVtcell{k} = VmaxweekavgCVt;
ninweekCVtcell{k} = ninweekCVt;
wknumCVtcell{k} = wknumCVt;
VmaxallcelltCV{k} = VmaxallCVt;
    
doseCVt = testing_vddata(:,2); % dose identifier in testing set
varCVt = testing_vddata (:,3); % viability identifier in testing set
cohort_numberCVt = testing_vddata(:,4); % cohort number in testing set (should only be one cohort)
nCVt = length(doseCVt); % finds number of data points total in testing set

options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
%define bounds for all parameters
            % m1 cen1 m2 cen 2 wk1sens wk1res etc...
paramslb = zeros([24 1]);
paramsub = Inf([24 1]);


params0 = [ .1; .1; .1; .1; .1; .1; .1; .1; .1; .1; .1; .1; 30; 30; 30; 30; 30; 30; 30; 30; 30; 30; 30; 30];
% outputs betaLD50 with first 12 paramters slope values, second 12 are
% center values in week order
[betaLD50t, resnormLD50t, residualsLD50st] = lsqnonlin(@fit_simpLD50,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    doseCVt,...
    varCVt,...
    wknumCVt,...
    VmaxallCVt);

betaLD50CVtesting(k,:) = betaLD50t;

end
%%
removedgaps = betaLD50CVtesting;
ind1 = removedgaps == 0.1;
ind2 = removedgaps == 30;
removedgaps(ind1) = NaN;
removedgaps(ind2) = NaN;
time = 1:1:12;
time = time';
for k = 1:num_cohorts
    cohort_params(:,1,k) = time;
    cohort_params(:,2,k)= removedgaps(k,1:12)';
    cohort_params(:,3,k) = removedgaps(k,13:24)';
    
end
%% Find naive LD50 of each cohort

for k = 1:num_cohorts
    paramslbn = zeros([2 1]);
    paramsubn = Inf([2 1]);
cohort_ind = naive_vddata(:,4) == cohorts(k);
cohort_naive= naive_vddata(cohort_ind,:);
dosenc = cohort_naive(:,2);
varnc = cohort_naive(:,3);

dose0ind = cohort_naive(:,2) == 0;
 Vmaxcohortnaivedata = cohort_naive(dose0ind,:);
 Vmaxnaiveavgcohort = mean(Vmaxcohortnaivedata(:,3));
params0n = [ .1; 30];
% outputs betaLD50 with first 12 paramters slope values, second 12 are
% center values in week order
[betaLD50naivecohort, resnormLD50naive, residualsLD50snaive] = lsqnonlin(@fit_simp1popnaive,...
    params0n,...
    paramslbn,...
    paramsubn,...
    options,...
    dosenc,...
    varnc,...
    Vmaxnaiveavgcohort);

betaLD50naivecohortall(k,:) = betaLD50naivecohort;
end


%% Plot Cohort Fits

for t = 1:12
for i = 1: num_cohorts
   all_cohorts(t,i) = cohort_params(t,3,i); % makes a matrix with each row a time and each column a cohort
 
end

end
for t = 1:12
 std_all_weeks(t) = nanstd(all_cohorts(t,:));
end

CI = 1.96.*std_all_weeks;

std_naive_LD50 = std(betaLD50naivecohortall(:,2));
CI_naive = 1.96*std_naive_LD50;

figure (3)
hold off
errorbar(time, betaLD50(13:24),CI,'--k', 'LineWidth', 1.5);
hold on
errorbar(0, betaLD50naiveall(2), CI_naive, '--k', 'LineWidth',1.5);
plot(time, cohort_params(:,3,1),'ro', 'LineWidth', 2);
plot(0, betaLD50naivecohortall(1,2), 'ro','LineWidth', 4);
plot( time, cohort_params(:,3,2),'mo', 'LineWidth', 2);
plot(0, betaLD50naivecohortall(2,2), 'mo','LineWidth',4);
plot( time, cohort_params(:,3,3),'yo', 'LineWidth', 2);
plot(0, betaLD50naivecohortall(3,2), 'yo','LineWidth',4);
plot(time, cohort_params(:,3,4),'go', 'LineWidth', 2);
plot(0, betaLD50naivecohortall(4,2), 'go','LineWidth',4);
plot(time, cohort_params(:,3,5),'co', 'LineWidth', 2);
plot(0, betaLD50naivecohortall(5,2), 'co','LineWidth',4);
plot(time, cohort_params(:,3,6),'bo', 'LineWidth', 2);
plot(0, betaLD50naivecohortall(6,2), 'bo','LineWidth',4);
xlabel('Time Post Treatment (weeks)', 'FontSize', 20);
xlim([-1 13])
ylim([0 100])
ylabel ('LD50 of each cohort', 'FontSize', 20);
title('Cohort Variation in LD50', 'FontSize', 20);

figure (4)
hold off
errorbar(time, betaLD50(13:24),CI,'--k', 'LineWidth', 1.5);
hold on
errorbar(0, betaLD50naiveall(2), CI_naive, '--k', 'LineWidth',1.5);
xlabel('Time Post Treatment (weeks)', 'FontSize', 20);
ylabel ('Best-Fit LD50', 'FontSize', 20);
xlim([-1 13])
ylim([0 100])
title('Parameter Estimates of LD50 with 95 % CI', 'FontSize', 20);




