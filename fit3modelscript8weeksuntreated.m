% This script takes data for all weeks and computes the average LD50s for
% each time point.  It plots the LD50s as a function of time. Then it does
% the same for each individual cohort, and plots the LD50s as a function of
% time on one graph

% Calculates Untreated values for only 8 weeks of data

clear all, close all, clc


vddata = load ('all_weeks_data_untreated.m'); % load concatenated data set (long Mx2 matrix
% all weeks added vertically, first column is week number
naive_vddata = load('naive_data.m');
%% Allocate data

week_ind = vddata(:,1) >8;
vddata = vddata(week_ind == 0, :);
coh_ind = vddata(:, 4) == 3;
vddata = vddata(coh_ind == 0, :);
%coh_ind8 = vddata(:,4) == 8;
%vddata = vddata(coh_ind8 == 0, :);
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
%% Extract cohort 8 from cell of testing sets

cohort_8 = cell_of_testingsets(4);
cohort_8 = cell2mat(cohort_8);

%% Find LD50 Naive

  options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
paramslbn = zeros([2 1]);
paramsubn = Inf([2 1]);

params0n = [ .1; 30];
% outputs betaLD50 with first 12 paramters slope values, second 12 are
% center values in week order
[betaLD50naive, resnormLD50naive, residualsLD50snaive] = lsqnonlin(@fit_simp1popnaive,...
    params0n,...
    paramslbn,...
    paramsubn,...
    options,...
    dosenc,...
    varnc,...
    Vmaxnaiveavg);
%% Bulk Population model constant over time

% Comes up with a model that has a single LD50 and slope for all the weeks
% data and allows only Vmax to vary for each week's data

initials1new = [ .05; 30]; % sets up initial values of beta1new

% call to lsqnonlin, which calls fit_simp1popnew which contains a
% difference vector in which the Vmax is replaced by different variables
% corresponding to the week. Here this parameter should be 1 every time
% (but its not)

[beta1unt_8wks, resnorm1, residuals1] = lsqnonlin(@fit_simp1popunw, initials1new,[0; 0],[ Inf; Inf],[],dose,var, wknum, Vmaxall);
X = (0:.1:250);
X = X'
v_model1=model1pop( dose, beta1unt_8wks, Vmaxweekavg );
v_model1plot = model1pop(X, beta1unt_8wks, Vmaxweekavg);
actualresiduals1 = v_model1-var;
sigma1 = std(residuals1)
error1 = sum(abs(residuals1))./n
chi_squared1 = sum((residuals1.^2)./v_model1)
RSS1 = sum(residuals1.^2);
p1 = 2;
AIC1 = n*log(RSS1./n) + 2*p1
BIC1 = -n*log(RSS1./n) + p1*log(n)
DOF1 = n-p1;
MSE1 = chi_squared1./ DOF1

plot(dose, var, '.', X, v_model1plot)

%% Perform Boot Strapping Estimates of Error in Parameters
[lowerlim1boot, upperlim1boot] = BSerrorinparams1(residuals1, dose, beta1unt_8wks, Vmaxweekavg, vddata);
[lowerlim1bootnaive, upperlim1bootnaive] = BSerrorinparams1naive(residualsLD50snaive, dosenc, betaLD50naive, naive_vddata, Vmaxnaiveavg);
%% Print Boot strapping limits
lowerlim1boot = lowerlim1boot
upperlim1boot = upperlim1boot
lowerlim1bootnaive = lowerlim1bootnaive
upperlim1bootnaive = upperlim1bootnaive
errorbarlength1boot = upperlim1boot-lowerlim1boot
errorbarlength1bootnaive = upperlim1bootnaive-lowerlim1bootnaive

%% Plot Bulk LD50
figure(1) % graphs fit to all weeks in red
hold off
plot (dose, var, 'o')
hold on
plot (X, v_model1plot, '-r', 'LineWidth', 3);
xlim([0,250]);
xlabel('dose (uM)','FontSize',18)
ylabel('viability','FontSize',18)
title ('Single Population Model All Data','FontSize',18)
set(gca,'LineWidth',1.5,'FontSize',18)

%% Dynamic Single Population Model: 8 weeks of data

options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
%define bounds for all parameters
            % m1 cen1 m2 cen 2 wk1sens wk1res etc...
paramslb = zeros([16 1]);
paramsub = Inf([16 1]);


params0 = [ .1; .1; .1; .1; .1; .1; .1; .1; 30; 30; 30; 30; 30; 30; 30; 30 ];
% outputs betaLD50 with first 12 paramters slope values, second 12 are
% center values in week order
[betaLD50unt_8wks, resnormLD50, residualsLD50] = lsqnonlin(@fit_simpLD508,...
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
%% Perform Boot Strapping Estimates of Error in Parameters Single Dynamic Model
[lowerlimLD50boot, upperlimLD50boot] = BSerrorinparamsLD50(residualsLD50, vddata, betaLD50unt_8wks);

%% Print Boot strapping limits
lowerlimLD50boot = lowerlimLD50boot
upperlimLD50boot = upperlimLD50boot
errorbarlengthLD50bootuntreated_8wks = upperlimLD50boot-lowerlimLD50boot


%% Plot of LD50s over time

time = (1:1:15);
errorbarsLD50 = errorbarlengthLD50bootuntreated';
figure (2)
% Plot of LD50s over time
%plot(time, betaLD50(16:30), 'o', 'LineWidth', 4)
errorbar(time, betaLD50unt(16:30), errorbarsLD50(16:30)/2, '-bo', 'LineWidth', 1.5)
hold on
errorbar(0, betaLD50naive(2), errorbarlength1bootnaive(2),'bo', 'LineWidth',1.5)
hold on
%plot(0, betaLD50naive(2), 'go', 'LineWidth',4)
xlabel ('Time Post Treatment', 'FontSize', 16)
ylabel ('Single Population LD50', 'FontSize', 16)
ylim ([20,80])
title('Single Dynamic Population Model for All Data','FontSize', 16)
set(gca,'LineWidth',1.5,'FontSize',20)

%% Estimate Monte Carlo Errors: Note needs to be rerun

upperlim = load('upperlim2.mat');
upperlim = struct2cell(upperlim);
upperlim = cell2mat(upperlim);
lowerlim = load('lowerlim2.mat');
lowerlim = struct2cell(lowerlim);
lowerlim = cell2mat(lowerlim);

CI_MC = upperlim -lowerlim;

avgMCfracts = mean(CI_MC(5:16));

CI_MC2 = vertcat(CI_MC(5:16), avgMCfracts-.003, avgMCfracts+.004, avgMCfracts+.005);


%% Model Fit Single Population Shifting in Time

bulk_model_LD50 = model_bulk_LD50( dose, betaLD50unt_8wks, Vmaxall, wk );
%bulk_model_LD50_plot = model_bulk_LD50_plot( X, betaLD50, Vmaxbyweek );
%% MC Error Estimation
[ lowerlimLD50, upperlimLD50, errorbarlengthLD50, betaLD50sim ] = MCerrorinparamsLD50( dose, bulk_model_LD50, betaLD50, kstd, wknum, Vmaxall);
CI_MC_LD50 = errorbarlengthLD50(16:30)./2;
%%

figure (3)
% Plot of LD50s over time
%plot(time, betaLD50(16:30), 'o', 'LineWidth', 4)
errorbar(time, betaLD50unt(16:30), errorbarsLD50(16:30)./2, 'b-', 'LineWidth', 2)
hold on
plot(0, betaLD50naive(2), 'bo', 'LineWidth',4)
xlabel ('Time Post Treatment', 'FontSize', 20)
ylabel ('Bulk Population LD50', 'FontSize', 20)
ylim ([20,80])
title('Bulk Population Dynamics for All Data','FontSize', 20)

figure(4)
hold off
% Plot of bulk model fits versus experiment
plot (dose, bulk_model_LD50,'o')
xlabel(' Dose Doxorubicin (uM)', 'FontSize', 20)
ylabel('Model Prediction of Viability', 'FontSize',20)

figure(5)
hold off
plot3(dose, vddata(:,1), bulk_model_LD50, 'ko', 'LineWidth',5)
hold on
plot3(dose, vddata(:,1), var, 'bo')
xlabel(' Dose Doxorubicin (uM)', 'FontSize', 12)
zlabel('Model Prediction and Experiment Viability', 'FontSize',12)
ylabel('Time (Weeks Post Treatment)','FontSize', 12)


%% Compute Statistics: Comparing Stagnant One Population to Dynamic One Population

actualresidualsLD50 = bulk_model_LD50-var;
sigmaLD50 = std(residualsLD50)
errorLD50 = sum(abs(residualsLD50))./n
chi_squaredLD50 = sum((residualsLD50.^2)./bulk_model_LD50)
RSSLD50 = sum(residualsLD50.^2);
pLD50 = 16;
AICLD50 = n*log(RSSLD50./n) + 2*pLD50
BICLD50 = -n*log(RSSLD50./n) + pLD50*log(n)
DOFLD50 = n-pLD50;
MSELD50 = chi_squaredLD50./ DOFLD50

% Comparison of one population changing in time with constant 1 population model
F_stat1vLD50alt = ((RSS1-RSSLD50)./(DOF1-DOFLD50))./(RSSLD50./(DOFLD50))
p1vLD50alt = 1-fcdf(F_stat1vLD50alt, DOF1-DOFLD50, DOFLD50)
%% Now find LD50s of the individual cohorts over time
% Outputs a 3D matrix containing the slope and LD50 fits over time for each
% cohort, will want to do the exact same for the fit to fsens to compute CI
% based on spread in that
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
cell_of_naive_testing_sets{k} = cohort_naive;

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

%% Extract naive data cohort 8

cohort_8_naive = cell_of_naive_testing_sets(4);
cohort_8_naive = cell2mat(cohort_8_naive);

%% Plot Cohort Fits for Dynamic Single Population Model

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

%%
for i= 1:6
cohort_params16_1_15(:,i) = cohort_params15(:,3,i)
end

time0s = betaLD50naivecohortall(:,2)';

cohort_params16 = vertcat( time0s, cohort_params16_1_15);

full_time = 0:1:15;
%%
figure (6)

betaLD50allunt = vertcat(betaLD50naive(2),betaLD50unt(16:30));
CI_untreated = horzcat(CI_naive, CI);
hold off
errorbar(full_time, betaLD50allunt,CI_untreated,'--k', 'LineWidth', 1.5);
hold on
plot(full_time, cohort_params16(:,1),'ro', 'LineWidth', 2);
plot(full_time, cohort_params16(:,2),'mo', 'LineWidth', 2);
plot( full_time, cohort_params16(:,3),'yo', 'LineWidth', 2);
plot(full_time, cohort_params16(:,4),'go', 'LineWidth', 2);
plot(full_time, cohort_params16(:,5),'co', 'LineWidth', 2);
plot(full_time, cohort_params16(:,6),'bo', 'LineWidth', 4);
xlabel('Time Post Treatment (weeks)', 'FontSize', 20);
xlim([-1 16])
ylim([0 100])
ylabel ('LD50 of each cohort', 'FontSize', 20);
title('Naive Cohort Variation in LD50', 'FontSize', 20);
legend('avg','3','4','6','8','9','10')

figure (5)
hold off
errorbar(time, betaLD50unt(16:30),CI,'--k', 'LineWidth', 1.5);
hold on
errorbar(0, betaLD50naive(2), CI_naive, 'ko', 'LineWidth',4);
xlabel('Time Post Treatment (weeks)', 'FontSize', 20);
ylabel ('Best-Fit LD50', 'FontSize', 20);
xlim([-1 16])
ylim([0 100])
title('Parameter Estimates of LD50 with 95 % CI', 'FontSize', 20);

% what a mess... compute confidence intervals

%% Two Population Model fractions shift in time
coeffs2 = load('coeffs2.mat')
coeffs2 = struct2cell(coeffs2);
coeffs2 = cell2mat(coeffs2);
options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
%define bounds for all parameters
            % m1 cen1 m2 cen 2 wk1sens wk1res etc...
paramslb = zeros([8 1]);
paramsub = [ 1; 1; 1; 1; 1; 1; 1; 1];

% sets up initial values of beta2new
params0 = [ .5; .5; .5; .5; .5; .5; .5; .5];
%  res2sq=ones([n 1]);
[beta2unt, resnorm2, residuals2] = lsqnonlin(@fit_simp2popunw8_untreated,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    dose,...
    var,...
    wknum,...
    Vmaxall,...
    coeffs2);

%% Compute and Output Two Population Statistics and Model Comparison
beta2_comb = vertcat(coeffs2, beta2unt);
v_model2allweeksnormed = model2popallweeksnormed( dose, beta2_comb, Vmaxbyweek, nsize);
actualresiduals2 = var- v_model2allweeksnormed;
sigma2 = std(residuals2)
error2 = sum(abs(residuals2))./n
chi_squared2 = sum((residuals2.^2)./v_model2allweeksnormed)
RSS2 = sum(residuals2.^2);
p2 = 12;
AIC2 = n*log(RSS2./n) + 2*p2
BIC2 = -n*log(RSS2./n) + p2*log(n)

DOF2 = n-p2;
MSE2 = chi_squared2./(DOF2)



figure(11)
hold off
plot3(dose, vddata(:,1), v_model2allweeksnormed, 'ro', 'LineWidth',5)
hold on
plot3(dose, vddata(:,1), var, 'bo')
xlabel(' Dose Doxorubicin (uM)', 'FontSize', 12)
zlabel('Model Prediction and Experiment Viability', 'FontSize',12)
ylabel('Time (Weeks Post Treatment)','FontSize', 12)


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

%% Find Naive Overall Fractions
; % loads 3 column vector of dose, viability, and cohort number
 
% simply want to average all data to find initial starting conditions of
% all cohorts
params0 = .5;
paramslb = 0;
paramsub = 1;
coeffs2 = beta2newall(1:4);
[beta2naive, resnorm2naive, residuals2naive] = lsqnonlin(@fit_simp2popnaiveunw,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    dosenc,...
    varnc,...
    coeffs2);

fsens0 = beta2naive;
fres0 = 1-beta2naive;
%% Find Naive model
v_model2allweeksnormednaive = model2popallweeksnormednaive( dose, coeffs2, beta2naive, Vmaxnaiveavg)

%% Perform Boot Strapping Estimates of Error in Parameters two population model
% Need to make this function and its naive counterpart
beta2_comb = vertcat(coeffs2, beta2unt);
[lowerlim2boot, upperlim2boot] = BSerrorinparams2(residuals2, beta2_comb, vddata);
%%
[lowerlim2bootnaive, upperlim2bootnaive] = BSerrorinparams2naive(residuals2naive, dosenc, beta2naive, naive_vddata, beta2newall, Vmaxnaiveavg);
%% Print Boot strapping limits
lowerlim2boot = lowerlim2boot
upperlim2boot = upperlim2boot
errorbarlength2bootuntreated = upperlim2boot-lowerlim2boot
lowerlim2bootnaive = lowerlim2bootnaive
upperlim2bootnaive = upperlim2bootnaive
errorbarlength2bootnaive = upperlim2bootnaive - lowerlim2bootnaive

%%
fres2 = zeros([15,1]);
fsens2 = beta2newall(5:19);

for i = 1:15
    fres2(i) = 1-fsens2(i);
end
time = zeros([1 15]);
time(1, :) = 1:15;

% Plot of the resistant and sensitive fractions
figure(6)
hold off
errorbar(time, fres2, errorbarlength2boot(5:19)./2, 'bo-', 'LineWidth', 1)
hold on
errorbar(0, fres0, errorbarlength2bootnaive./2, 'bo', 'LineWidth', 1)
xlim([ 0 15])
ylim([ 0 1])
xlabel('Time (Weeks Post Treatment)', 'FontSize', 16)
ylabel('Fraction of Cells', 'FontSize',16)
title ('Fraction of Resistant Cells Two Population Model', 'FontSize', 14)
hold off

%% Model Fit for 3 example weeks on one graph
v_model2_3 = model2pop16p( X, beta2newall, Vmaxbyweek, 3);
v_model2_6 = model2pop16p( X, beta2newall, Vmaxbyweek, 6);
v_model2_10 = model2pop16p( X, beta2newall, Vmaxbyweek, 10);
time(1, :) = 1:15;

figure(7) % looks at example weeks 3, 6, and 10 WPT against all experimental data
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

%% Run exact same cohort analysis on two population model

% Find Fres and Fsens of the individual cohorts over time
% Outputs a 3D matrix containing the slope and fsens fits over time for each
% cohort, will want to do the exact same for the fit to fsens to compute CI
% based on spread in that
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
  paramslb = zeros([19 1]);
paramsub = [ Inf; Inf; Inf; Inf; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1;];

% sets up initial values of beta2new
params0 = [ .1; 17; .04; 50; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5];
%  res2sq=ones([n 1]);
[beta2newt, resnorm2, residuals2] = lsqnonlin(@fit_simp2popunw15,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    doseCVt,...
    varCVt,...
    wknumCVt,...
    VmaxallCVt);             
               




beta2newCVtesting15(k,:) = beta2newt;

end
%%
removedgaps2 = beta2newCVtesting15;
ind1 = removedgaps2 == 0.5;
removedgaps2(ind1) = NaN;
removedgaps2= removedgaps2';
%%
time = 1:1:15;
time = time';
for k = 1:num_cohorts
    cohort_params2_15(:,1,k) = time;
    cohort_params2_15(:,2,k)= removedgaps2(5:19,k); 
    coeffs2cohort(:,k) = removedgaps2(1:4,k);
end

%% Find naive fsens of each cohort

for k = 1:num_cohorts
    paramslbn = 0;
    paramsubn = 1;
cohort_ind = naive_vddata(:,4) == cohorts(k);
cohort_naive= naive_vddata(cohort_ind,:);
dosenc = cohort_naive(:,2);
varnc = cohort_naive(:,3);

dose0ind = cohort_naive(:,2) == 0;
 Vmaxcohortnaivedata = cohort_naive(dose0ind,:);
 Vmaxnaiveavgcohort = mean(Vmaxcohortnaivedata(:,3));
 coeffs2nc = coeffs2cohort(:,k);
params0n = .6;
% outputs betaLD50 with first 12 paramters slope values, second 12 are
% center values in week order
[beta2naivecohort, resnorm2nc, residuals2nc] = lsqnonlin(@fit_simp2popnaiveunw,...
    params0n,...
    paramslbn,...
    paramsubn,...
    options,...
    dosenc,...
    varnc,...
    coeffs2nc);

beta2naivecohortall(k,:) = beta2naivecohort;
end


%% Plot Cohort Fits for Dynamic Single Population Model

for t = 1:15
for i = 1: num_cohorts
   all_cohorts2(t,i) = cohort_params2_15(t,2,i); % makes a matrix with each row a time and each column a cohort
 
end

end
for t = 1:15
 std_all_weeks2(t) = nanstd(all_cohorts2(t,:)); % row with each column containing std by week
end
for b = 1:4
std_params(b) = std(coeffs2cohort(b,:))
end

CI2 = 1.96.*std_all_weeks2;
CIparams = 1.96.*std_params

std_naive_2 = std(beta2naivecohortall);
CI_naive_2 = 1.96*std_naive_2;
%%
figure (8)
hold off
errorbar(time, fres2, CI2,'--k', 'LineWidth', 1.5);
hold on
errorbar(0,1- beta2naive, CI_naive_2, 'ko', 'LineWidth',4);
plot(time, 1-cohort_params2_15(:,2,1),'ro', 'LineWidth', 2);
plot(0, 1-beta2naivecohortall(1), 'ro','LineWidth', 4);
plot( time, 1-cohort_params2_15(:,2,2),'mo', 'LineWidth', 2);
plot(0, 1-beta2naivecohortall(2), 'mo','LineWidth', 4);
plot( time, 1-cohort_params2_15(:,2,3),'yo', 'LineWidth', 2);
plot(0, 1-beta2naivecohortall(3), 'yo','LineWidth', 4);
plot(time, 1-cohort_params2_15(:,2,4),'go', 'LineWidth', 2);
plot(0, 1-beta2naivecohortall(4), 'go','LineWidth', 4);
plot(time, 1-cohort_params2_15(:,2,5),'co', 'LineWidth', 2);
plot(0, 1-beta2naivecohortall(5), 'co','LineWidth', 4);
plot(time, 1-cohort_params2_15(:,2,6),'bo', 'LineWidth', 2);
plot(0, 1-beta2naivecohortall(6), 'bo','LineWidth', 4);
xlabel('Time Post Treatment (weeks)', 'FontSize', 20);
xlim([-1 16])
ylim([0 1])
ylabel ('Fraction Resistant each cohort', 'FontSize', 20);
title('Cohort Variation in Resistant Fraction', 'FontSize', 20);
%%
figure (9)
hold off
errorbar(time, 1- beta2newall(4:19),CI2,'--k', 'LineWidth', 1.5);
hold on
errorbar(0, beta2naive(2), CI_naive, 'ko', 'LineWidth',4);
xlabel('Time Post Treatment (weeks)', 'FontSize', 20);
ylabel ('Best-Fit LD50', 'FontSize', 20);
xlim([-1 16])
ylim([0 100])
title('Parameter Estimates of Fraction Resistant with 95 % CI', 'FontSize', 20);

% what a mess... compute confidence intervals


%% Two Population Model Testing and Training Models

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
% does the same but across all doses given by X
end

coeffs2train= beta2newCVtraining(1:num_cohorts, 1:4);
% k by 4 matrix of slope and center paramters from each training set model
% for clarity from fractional estimate parameters

%% 
%% Two Population General Model Testing Sets
% Now compare to the fit of a single cohort given by the testing set
% For each k, extract testing set 
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
nCVt = length(dose); % finds number of data points total in testing set

options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
%define bounds for all parameters
            % m1 cen1 m2 cen 2 wk1sens wk1res etc...
paramslb = zeros([15 1]);
paramsub = ones([15 1]);

% sets up initial values of beta2new
params0 = [ .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; ,5; .5; .5];
%  res2sq=ones([n 1]);

[beta2new15p, resnorm2, residuals2] = lsqnonlin(@fit_simp2pop15pnormfix,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    doseCVt,...
    varCVt,...
    wknumCVt,...
    VmaxallCVt,...
    coeffs2train(k,:));
%coeffs2(k,:))
beta2newCVtesting(k,:) = beta2new15p; %% This creates the testing predictions of fsens

end

removedgapsb2newt = beta2newCVtesting;
timeind = removedgapsb2newt == 0.5;
removedgapsb2newt(timeind) = 0;
time = 1:1:15;
time = time';
for k = 1:num_cohorts
    fractionvt(:,1,k) = time;
    fractionvt(:,2,k)= removedgapsb2newt(k,:)';
end
fractionvt(:,3,:) = 1- fractionvt(:,2,:);
% this creates a 3d vector where in each k dimension is the fraction fits
% of the testing data to the beta2avg parameters

%% Compare and Compute testing set estimates of resistant and sensitive fractions

for k = 1:num_cohorts
    fractsmanip = fractionvt(:,2,k);
%     fractsmanip = fractsmanip';

    updatedfractions = fractionvt(:,:,k);
   
ind = updatedfractions(:,2) == 0;

updatedfractions(ind,:) = [];

celloffractionvtime{k} = updatedfractions;


testfractsens = updatedfractions(:,2);
sensfractsavgmanip = sensfractsavg;
sensfractsavgmanip(ind) = [];

celloferrorinfracts{k} = testfractsens-sensfractsavgmanip;


errorinfracts = testfractsens-sensfractsavgmanip;

avgerror(k) = mean(abs(errorinfracts));
end


%% Plots to compare testing set to 95 % CI of training sets


fractsvtime1 = celloffractionvtime(1);
fractsvtime1 = cell2mat(fractsvtime1);

fractsvtime2 = celloffractionvtime(2);
fractsvtime2 = cell2mat(fractsvtime2);

fractsvtime3 = celloffractionvtime(3);
fractsvtime3 = cell2mat(fractsvtime3);

fractsvtime4 = celloffractionvtime(4);
fractsvtime4 = cell2mat(fractsvtime4);

fractsvtime5 = celloffractionvtime(5);
fractsvtime5 = cell2mat(fractsvtime5);

fractsvtime6 = celloffractionvtime(6);
fractsvtime6 = cell2mat(fractsvtime6);
% plot all testing set cohorts against final model

figure(6)
hold off
plot(time, resfractsavg, '--k', 'LineWidth', 4)
%errorbar(time, fres2, errorbarlength2(5:16)./2, 'b', 'LineWidth', 3)
hold on
plot(1:12, resfracts95CI(:,1),'-r', 'LineWidth', 4)
plot(1:12, resfracts95CI(:,2),'-r', 'LineWidth', 4)
xlim([ 1 12])
ylim( [0 1])
plot(fractsvtime1(:,1), fractsvtime1(:,3), '-r', 'LineWidth', 2)
plot( 0, 1-0.8050, 'go', 'LineWidth', 4)
plot(fractsvtime2(:,1), fractsvtime2(:,3), '-m', 'LineWidth', 2)
plot( 0, 1-0.80840, 'co', 'LineWidth', 4)
plot(fractsvtime3(:,1), fractsvtime3(:,3), '-y','LineWidth', 2)
plot( 0, 1-0.3061, 'bo', 'LineWidth', 4)
plot(fractsvtime4(:,1), fractsvtime4(:,3), '-g', 'LineWidth', 2)
plot(fractsvtime5(:,1), fractsvtime5(:,3), '-c', 'LineWidth', 2)
plot(fractsvtime6(:,1), fractsvtime6(:,3), '-b','LineWidth', 2)
%errorbar(time, fsens2, errorbarlength2(5:16)./2, 'g', 'LineWidth',3)
xlabel('Weeks Post Treatment', 'FontSize', 24)
ylabel('Fraction of Cells', 'FontSize',24)
title('Cohort Variation in Resistant Fraction', 'FontSize', 18)

figure(7)
hold off
plot(time, sensfractsavg, '--k', 'LineWidth',4)
hold on
plot(1:12, beta295CI(5:16, 1),'-r', 'LineWidth', 4)
plot(1:12, beta295CI(5:16, 2),'-r', 'LineWidth', 4)
plot(fractsvtime1(:,1), fractsvtime1(:,2), '-r', 'LineWidth', 2)
plot( 0, 0.8050, 'go', 'LineWidth', 4)
plot(fractsvtime2(:,1), fractsvtime2(:,2), '-m', 'LineWidth', 2)
plot( 0, 0.8084, 'co', 'LineWidth', 4)
plot(fractsvtime3(:,1), fractsvtime3(:,2), '-y','LineWidth', 2)
plot( 0, 0.3061, 'bo', 'LineWidth', 4)
xlim([1 12])
ylim([0 1])
plot(fractsvtime4(:,1), fractsvtime4(:,2), '-g', 'LineWidth', 2)
plot(fractsvtime5(:,1), fractsvtime5(:,2), '-c', 'LineWidth', 2)
plot(fractsvtime6(:,1), fractsvtime6(:,2), '-b','LineWidth', 2)
%errorbar(time, fsens2, errorbarlength2(5:16)./2, 'g', 'LineWidth',3)
xlabel('Weeks Post Treatment', 'FontSize', 24)
ylabel('Fraction of Cells', 'FontSize',24)
title('Cohort Variation in Sensitive Fraction', 'FontSize', 18)


%This returns a vector with each cohort in the z direction and each column
%is wpt, fsens, fres and each row is the value corresponding to that column




