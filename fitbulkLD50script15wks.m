% This script takes data for all weeks and computes the average LD50s for
% each time point.  It plots the LD50s as a function of time. Then it does
% the same for each individual cohort, and plots the LD50s as a function of
% time on one graph

clear all, close all, clc


vddata = load ('allweeksdata12_13_15.m'); % load concatenated data set (long Mx2 matrix
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

%% Bulk Population model constant over time

% Comes up with a model that has a single LD50 and slope for all the weeks
% data and allows only Vmax to vary for each week's data

initials1new = [ .05; 30]; % sets up initial values of beta1new

% call to lsqnonlin, which calls fit_simp1popnew which contains a
% difference vector in which the Vmax is replaced by different variables
% corresponding to the week. Here this parameter should be 1 every time
% (but its not)

[beta1new, resnorm1, residuals1] = lsqnonlin(@fit_simp1popunw, initials1new,[0; 0],[ Inf; Inf],[],dose,var, wknum, Vmaxall);
X = (0:.1:250);
v_model1=model1pop( dose, beta1new, Vmaxweekavg );
v_model1plot = model1pop(X, beta1new, Vmaxweekavg);
actualresiduals1 = v_model1-var;
%chi2=sum((residual.^2)./S_t);
sigma1 = std(residuals1)
error1 = sum(abs(residuals1))./n
chi_squared1 = sum((residuals1.^2)./v_model1)
RSS1 = sum(residuals1.^2);
p1 = 2;
AIC1 = n*log(RSS1./n) + 2*p1
BIC1 = -n*log(RSS1./n) + p1*log(n)
DOF1 = n-p1;
MSE1 = chi_squared1./ DOF1

%% Graph of Bulk LD50
figure(1) % graphs fit to all weeks in red
hold off
plot (dose, var, 'o')
hold on
plot (X, v_model1plot, '-r', 'LineWidth', 3);
xlim([0,250]);
xlabel('dose (uM)','FontSize',18)
ylabel('viability','FontSize',18)
title ('Single Population Model','FontSize',18)

%% 15 weeks of data: one population model shifting in time

options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
%define bounds for all parameters
            % m1 cen1 m2 cen 2 wk1sens wk1res etc...
paramslb = zeros([30 1]);
paramsub = Inf([30 1]);


params0 = [ .1; .1; .1; .1; .1; .1; .1; .1; .1; .1; .1; .1; .1; .1; .1; 30; 30; 30; 30; 30; 30; 30; 30; 30; 30; 30; 30; 30; 30; 30];
% outputs betaLD50 with first 12 paramters slope values, second 12 are
% center values in week order
[betaLD50, resnormLD50, residualsLD50] = lsqnonlin(@fit_simpLD5015,...
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
[betaLD50naive, resnormLD50naive, residualsLD50snaive] = lsqnonlin(@fit_simp1popnaive,...
    params0n,...
    paramslbn,...
    paramsubn,...
    options,...
    dosen,...
    varn,...
    Vmaxnaiveavg);


%% Plot of LD50s

time(1, :) = 1:15;
figure (2)
% Plot of LD50s over time
plot(time, betaLD50(16:30), 'o', 'LineWidth', 4)
hold on
plot(0, betaLD50naive(2), 'go', 'LineWidth',4)
xlabel ('Time Post Treatment', 'FontSize', 20)
ylabel ('Bulk Population LD50', 'FontSize', 20)
ylim ([20,80])
title('Bulk Population Dynamics for All Data','FontSize', 20)

%% Model Fit

bulk_model_LD50 = model_bulk_LD50( dose, betaLD50, Vmaxall, wk );

figure(3)
% Plot of bulk model fits versus experiment
plot (dose, bulk_model_LD50,'o')
xlabel(' Dose Doxorubicin (uM)', 'FontSize', 20)
ylabel('Model Prediction of Viability', 'FontSize',20)

%% Compute All Statistics

actualresidualsLD50 = bulk_model_LD50-var;
sigmaLD50 = std(residualsLD50)
errorLD50 = sum(abs(residualsLD50))./n
chi_squaredLD50 = sum((residualsLD50.^2)./bulk_model_LD50)
RSSLD50 = sum(residualsLD50.^2);
pLD50 = 30;
AICLD50 = n*log(RSSLD50./n) + 2*pLD50
BICLD50 = -n*log(RSSLD50./n) + pLD50*log(n)
DOFLD50 = n-pLD50;
MSELD50 = chi_squaredLD50./ DOFLD50

% Comparison of one population changing in time with constant 1 population model
F_stat1vLD50alt = ((RSS1-RSSLD50)./(DOF1-DOFLD50))./(RSSLD50./(DOFLD50))
p1vLD50alt = 1-fcdf(F_stat1vLD50alt, DOF1-DOFLD50, DOFLD50)

%% Two Population Model fractions shift in time

options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
%define bounds for all parameters
            % m1 cen1 m2 cen 2 wk1sens wk1res etc...
paramslb = zeros([19 1]);
paramsub = [ Inf; Inf; Inf; Inf; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1;];

% sets up initial values of beta2new
params0 = [ .1; 17; .04; 50; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5];
%  res2sq=ones([n 1]);
[beta2newall, resnorm2, residuals2] = lsqnonlin(@fit_simp2popunw15,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    dose,...
    var,...
    wknum,...
    Vmaxall);

%% Compute and Output Two Population Statistics and Model Comparison

v_model2allweeksnormed = model2popallweeksnormed( dose, beta2newall, Vmaxbyweek, nsize);
actualresiduals2 = var- v_model2allweeksnormed;
sigma2 = std(residuals2)
error2 = sum(abs(residuals2))./n
chi_squared2 = sum((residuals2.^2)./v_model2allweeksnormed)
RSS2 = sum(residuals2.^2);
p2 = 19;
AIC2 = n*log(RSS2./n) + 2*p2
BIC2 = -n*log(RSS2./n) + p2*log(n)

DOF2 = n-19;
MSE2 = chi_squared2./(DOF2)

% Comparison of two population model with bulk constant time model
F_stat1v2alt = ((RSS1-RSS2)./(DOF1-DOF2))./(RSS2./(DOF2))
p1v2alt = 1-fcdf(F_stat1v2alt, DOF1-DOF2, DOF2)

% Comparison of one population changing in time with two population model
F_stat2vLD50alt = ((RSS2-RSSLD50)./(DOF2-DOFLD50))./(RSSLD50./(DOFLD50))
p2vLD50alt = 1-fcdf(F_stat2vLD50alt, DOF2-DOFLD50, DOFLD50)
% Here 2 population model is the simpler model
%If F~1 simpler model is adequate
% If F>1 the more complex model is better, or random error led to a better
% fit with the complex model
%%
fres2 = zeros([15,1]);
fsens2 = beta2newall(5:19);

for i = 1:15
    fres2(i) = 1-fsens2(i);
end
time = zeros([1 15]);
time(1, :) = 1:15;

% Find error bars
% Note that upper and lower limits of parameter valyes found using Monte
% Carlo approximation for error in parameters based on standard deviation
% in viability in technical replicates (kstd function). Code not included
% here as is very time consuming.

lowerlim2 = load('lowerlim2.mat');
lowerlim2 = struct2cell(lowerlim2);
lowerlim2 = cell2mat(lowerlim2);
upperlim2 = load('upperlim2.mat');
upperlim2 = struct2cell(upperlim2);
upperlim2 = cell2mat(upperlim2);

errorbarlength2 = upperlim2-lowerlim2;

% Plot of the resistant and sensitive fractions
figure(4)
hold off
plot(time, fres2, 'b', 'LineWidth', 4)
hold on
%errorbar(time, fres2, errorbarlength2(5:16), 'b', 'LineWidth', 3)
plot(time, fsens2, 'g', 'LineWidth',4)
errorbar(time, fsens2, errorbarlength2(5:16), 'g', 'LineWidth',3)
xlim([ 1 12])
ylim([ 0 1])
xlabel('Weeks Post Treatment', 'FontSize', 24)
ylabel('Fraction of Cells', 'FontSize',24)
title ('Fraction of Resistant Cells and Sensitive Cells vs. Time Post Treatment from All Data', 'FontSize', 14)
hold off

%% Displays Fit for 3 example weeks on one graph
v_model2_3 = model2pop16p( X, beta2newall, Vmaxbyweek, 3);
v_model2_6 = model2pop16p( X, beta2newall, Vmaxbyweek, 6);
v_model2_10 = model2pop16p( X, beta2newall, Vmaxbyweek, 10);
time(1, :) = 1:15;

figure(3) % looks at example weeks 3, 6, and 10 WPT against all experimental data
hold off
plot (dose(sum(nsize(1:2))+1:sum(nsize(1:3))), var(sum(nsize(1:2))+1:sum(nsize(1:3))), 'ro')
hold on
plot (dose(sum(nsize(1:5))+1:sum(nsize(1:6))), var(sum(nsize(1:5))+1:sum(nsize(1:6))), 'go')
plot (dose(sum(nsize(1:9))+1:sum(nsize(1:10))), var(sum(nsize(1:9))+1:sum(nsize(1:10))), 'bo')
plot (X, v_model2_3, '-r', 'LineWidth', 3);
plot (X, v_model2_6, '-g', 'LineWidth', 3);
plot (X, v_model2_10, '-b', 'LineWidth', 3);
xlim ([0 250])
xlabel('dose (uM)','FontSize',24)
ylabel('viability','FontSize',24)
title ('Two Population Model 3, 6 & 10 Weeks Post Treatment','FontSize',18)

%% Create Training and Testing Sets to Input into all Three Models

% Identifies a cohort and extracts it to create k training sets with one
% cohort removed
for k=1:num_cohorts
    data_maniptrain = vddata; % sets training dating = vddata
    TF1 = data_maniptrain(:,4)==cohorts(k); % sets condition
    % removes those with that condition
    data_maniptrain(TF1,:) = [] ;
    cell_of_trainingsets{k} = data_maniptrain; % puts training set into cell array             
end   
 % Identifies cohort and sets tat cohort as the testing set (to be compared
 % to the training set with that cohort removed)
for k=1:num_cohorts
    data_maniptest = vddata;
    index1 = data_maniptest(:,4) == cohorts(k);
    testingset = data_maniptest(index1,:);
    cell_of_testingsets{k} = testingset;   %puts testing set into cell array          
end   
%% Two Population Model Testing and Training Models

X = X';% creates X values to be used to plot model parameters

% Training sets allocated into cells (because of all different sizes)
% Converts to matrix, finds all necessary Vmax data for each training set,
% and computes model parameters based on training set data
for k = 1:num_cohorts
    training_vddata = cell_of_trainingsets(k);
    training_vddata = cell2mat(training_vddata);
    
[ VmaxbyweekCV, VmaxweekavgCV, ninweekCV, wknumCV, VmaxallCV] = findVmaxandsize(training_vddata);
nsizeCV = wknumCV(:,2);

VmaxbyweekCVcell{k} = VmaxbyweekCV;
VmaxweekavgCVcell{k} = VmaxweekavgCV;
ninweekCVcell{k} = ninweekCV;
wknumCVcell{k} = wknumCV;
VmaxallcellCV{k} = VmaxallCV;
    
doseCV = training_vddata(:,2); % identifies dose in vddata of trainng set
varCV = training_vddata (:,3); % identifies viability in vddata of training set
cohort_numberCV = training_vddata(:,4); % identifies the cohort number in vddata of training set
nCV = length(doseCV); % finds number of data points total in training set

options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
%define bounds for all parameters
            % m1 cen1 m2 cen 2 wk1sens wk1res etc...
paramslb = zeros([19 1]);
paramsub = [ Inf; Inf; Inf; Inf; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1];

% sets up initial values of beta2new
params0 = [ .1; 17; .04; 50; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5];
%  res2sq=ones([n 1]);
[beta2new, resnorm2, residuals2] = lsqnonlin(@fit_simp2popunw15,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    doseCV,...
    varCV,...
    wknumCV,...
    VmaxallCV);

beta2newCVtraining(k,:) = beta2new;
 % makes k by 16 parameter sets based on each training model
 
v_model2trainingsetcell{k} =model2popallweeksnormed( doseCV, beta2new,VmaxbyweekCV, nsizeCV );
% uses each training set data and model paramters to create viability
% vector which represents that training models prediction
v_model2trainingsetplotcell{k} =model2popallweeksnormed( X, beta2new, VmaxbyweekCV,nsizeCV );
% does the same but across all doses given by X
end

coeffs2= beta2newCVtraining(1:num_cohorts, 1:4);
% k by 4 matrix of slope and center paramters from each training set model
% for clarity from fractional estimate parameters
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
    
 [ VmaxbyweekCVt, VmaxweekavgCVt, ninweekCVt, wknumCVt, VmaxallCVt] = findVmaxandsizetest15(testing_vddata);
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
paramslb = zeros([30 1]);
paramsub = Inf([30 1]);


params0 = [ .1; .1; .1; .1; .1;.1;.1; .1; .1; .1; .1; .1; .1; .1; .1; 30; 30; 30; 30; 30; 30; 30; 30; 30; 30; 30; 30; 30; 30; 30];
% outputs betaLD50 with first 12 paramters slope values, second 12 are
% center values in week order
[betaLD50t, resnormLD50t, residualsLD50st] = lsqnonlin(@fit_simpLD5015,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    doseCVt,...
    varCVt,...
    wknumCVt,...
    VmaxallCVt);

betaLD50CVtesting15(k,:) = betaLD50t;

end
%%
removedgaps = betaLD50CVtesting15;
ind1 = removedgaps == 0.1;
ind2 = removedgaps == 30;
removedgaps(ind1) = NaN;
removedgaps(ind2) = NaN;
time = 1:1:15;
time = time';
for k = 1:num_cohorts
    cohort_params15(:,1,k) = time;
    cohort_params15(:,2,k)= removedgaps(k,1:15)';
    cohort_params15(:,3,k) = removedgaps(k,16:30)';
    
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

for t = 1:15
for i = 1: num_cohorts
   all_cohorts(t,i) = cohort_params15(t,3,i); % makes a matrix with each row a time and each column a cohort
 
end

end
for t = 1:15
 std_all_weeks(t) = nanstd(all_cohorts(t,:));
end

CI = 1.96.*std_all_weeks;

std_naive_LD50 = std(betaLD50naivecohortall(:,2));
CI_naive = 1.96*std_naive_LD50;

figure (3)
hold off
errorbar(time, betaLD50(16:30),CI,'--k', 'LineWidth', 1.5);
hold on
errorbar(0, betaLD50naive(2), CI_naive, 'ko', 'LineWidth',4);
plot(time, cohort_params15(:,3,1),'ro', 'LineWidth', 2);
plot(0, betaLD50naivecohortall(1,2), 'ro','LineWidth', 4);
plot( time, cohort_params15(:,3,2),'mo', 'LineWidth', 2);
plot(0, betaLD50naivecohortall(2,2), 'mo','LineWidth', 4);
plot( time, cohort_params15(:,3,3),'yo', 'LineWidth', 2);
plot(0, betaLD50naivecohortall(3,2), 'yo','LineWidth', 4);
plot(time, cohort_params15(:,3,4),'go', 'LineWidth', 2);
plot(0, betaLD50naivecohortall(4,2), 'go','LineWidth', 4);
plot(time, cohort_params15(:,3,5),'co', 'LineWidth', 2);
plot(0, betaLD50naivecohortall(5,2), 'co','LineWidth', 4);
plot(time, cohort_params15(:,3,6),'bo', 'LineWidth', 2);
plot(0, betaLD50naivecohortall(6,2), 'bo','LineWidth', 4);
xlabel('Time Post Treatment (weeks)', 'FontSize', 20);
xlim([-1 16])
ylim([0 100])
ylabel ('LD50 of each cohort', 'FontSize', 20);
title('Cohort Variation in LD50', 'FontSize', 20);

figure (4)
hold off
errorbar(time, betaLD50(16:30),CI,'--k', 'LineWidth', 1.5);
hold on
errorbar(0, betaLD50naive(2), CI_naive, 'ko', 'LineWidth',4);
xlabel('Time Post Treatment (weeks)', 'FontSize', 20);
ylabel ('Best-Fit LD50', 'FontSize', 20);
xlim([-1 16])
ylim([0 100])
title('Parameter Estimates of LD50 with 95 % CI', 'FontSize', 20);

% what a mess... compute confidence intervals

%% Add in two Population Model 15 weeks




