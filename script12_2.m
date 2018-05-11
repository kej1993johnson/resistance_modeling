clear all, close all, clc
% This script takes data for all weeks, computes the one and two population models,
% performs the relevant statistics for model comparison, and removes cohorts for training and
% testing


vddata = load ('allweeksdata11_22.m'); % load concatenated data set (long Mx2 matrix
% all weeks added vertically, first column is week number
%%

dose = vddata (:,2);
var = vddata (:,3);
cohort_number = vddata(:,4);
n = length(dose); % finds number of data points total
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

%% Finds error based on fit and technical error to be used in weighting cost functions

[k1,k2] = fitabsoluteresiduals( dose, var, wknum, Vmaxall);
% This function fits all of the data to the two models, the finds the average residual
% value of each dose point (given in a range), fits these absolute residuals to a 
% model, and inverse of the absolute value of the residual is used to weight the
% cost function in our next fit

[kstd, Wstd, std_dev_model_fun] = fittechnicalerrors(dose);

% Fits the same exponential function except instead of using residuals uses
% standard deviation in technical replicates at each dose. Since doses are the same,
% binning as in above isn't required and function is much simpler. One
% issue is that this may not be extrapolated to be the proper error as a
% function of dose

%% One Population General Model
% Comes up with a model that has a single LD50 and slope for all the weeks
% data and allows only Vmax to vary for each week's data

initials1new = [ .05; 30]; % sets up initial values of beta1new

% call to lsqnonlin, which calls fit_simp1popnew which contains a
% difference vector in which the Vmax is replaced by different variables
% corresponding to the week. Here this parameter should be 1 every time
% (but its not)

[beta1new, resnorm1, residuals1] = lsqnonlin(@fit_simp1popabsresdens, initials1new,[0; 0],[ Inf; Inf],[],dose,var, k2, wknum, Vmaxall);

%%
% calculates AIC value based on general model, not sure about p value here
% since we're technically fitting this data to 14 betas (slope, cen, +12
% weeks Vmaxes)


disp (' Bulk slope and x-center')
new_coeffs_1pop = beta1new

%new_y_1pop = (Vmaxavg./( 1 + exp(beta1new(1).*(X - beta1new(2))))) ;
X = (0:.1:250);
v_model1=model1pop( dose, beta1new, Vmaxweekavg );
v_model1plot = model1pop(X, beta1new, Vmaxweekavg);
actualresiduals1 = v_model1-var;
%chi2=sum((residual.^2)./S_t);
sigma1 = std(actualresiduals1)
error1 = sum(abs(actualresiduals1))./n
chi_squared1 = sum((actualresiduals1.^2)./v_model1)
RSS1 = sum(actualresiduals1.^2);
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
xlim([0,200]);
xlabel('dose (uM)','FontSize',18)
ylabel('viability','FontSize',18)
title ('Single Population Model','FontSize',18)

%% Two Population General Model
% Creates a model with one set of slope and LD50s for a sensitive and
% resistant population, then fits each weeks data to a different fraction
% that corresponds to that LD50 and slope
% call to lsqnonlin, which calls fit_simp2popnew which contains a
% difference vector in which the f1(fraction of cells in first population
% is replaced by different variables corresponding to the week. Here this 
% parameter should be some fraction of the total that is sensitive for each
% week. Used two fractions and they don't necessarily add to 1 -> need to
% fix this going forward

% LSQNONLIN solves a particular class of minimization problems
% of the form min sum[(f(beta2new;dose_i) - var_i)^2] where the dose_i and
% var_i are data values and beta2new is the vector of unknowns. 
options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
%define bounds for all parameters
            % m1 cen1 m2 cen 2 wk1sens wk1res etc...
paramslb = zeros([16 1]);
paramsub = [ Inf; Inf; Inf; Inf; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1];

% sets up initial values of beta2new
params0 = [ .1; 17; .04; 50; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5];
%  res2sq=ones([n 1]);
[beta2newall, resnorm2, residuals2] = lsqnonlin(@fit_simp2popabsresdensnormed,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    dose,...
    var,...
    k2,...
    wknum,...
    Vmaxall);
%% Find MC Error in Parameters based on technical error function

[ lowerlimall, upperlimall, errorbarlengthall, beta2newsim ] = MCerrorinparams2( dose, beta2newall, Vmaxbyweek, k2, kstd, wknum, Vmaxall);
%% Compute and Output Two Population Statistics and Model Comparison

disp ('Sensitive and Resistant Population Slope and LD50s')
coeffs2all = beta2newall(1:4)
v_model2allweeksnormed = model2popallweeksnormed( dose, beta2newall, Vmaxbyweek, nsize);
actualresiduals2 = var- v_model2allweeksnormed;
sigma2 = std(actualresiduals2)
error2 = sum(abs(actualresiduals2))./n
chi_squared2 = sum((actualresiduals2.^2)./v_model2allweeksnormed)
RSS2 = sum(actualresiduals2.^2);
p2 = 16;
AIC2 = n*log(RSS2./n) + 2*p2
BIC2 = -n*log(RSS2./n) + p2*log(n)

coeffs2all = beta2newall(1:4)

DOF1 = n-2;
DOF2 = n-16;
MSE2 = chi_squared2./(DOF2)

F_stat1v2 = ((chi_squared1-chi_squared2)./(DOF1-DOF2))./(chi_squared2./(DOF2))
p1v2 = 1-fcdf(F_stat1v2, DOF1-DOF2, DOF2)
F_stat1v2alt = ((RSS1-RSS2)./(DOF1-DOF2))./(RSS2./(DOF2))
p1v2alt = 1-fcdf(F_stat1v2alt, DOF1-DOF2, DOF2)
%If F~1 simpler model is adequate
% If F>1 the more complex model is better, or random error led to a better
% fit with the complex model
%%
fres2 = zeros( [12,1]);
fsens2 = beta2newall(5:16);

for i = 1:12
    fres2(i) = 1-fsens2(i);
end
time = zeros([1 12]);
time(1, :) = 1:12;

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
figure(2)
hold off
plot(time, fres2, 'b', 'LineWidth', 4)
hold on
errorbar(time, fres2, errorbarlength2(5:16), 'b', 'LineWidth', 3)
plot(time, fsens2, 'g', 'LineWidth',4)
errorbar(time, fsens2, errorbarlength2(5:16), 'g', 'LineWidth',3)
xlim([ 1 12])
ylim([ 0 1])
xlabel('Weeks Post Treatment', 'FontSize', 24)
ylabel('Fraction of Cells', 'FontSize',24)
%title ('Fraction of Resistant Cells and Sensitive Cells vs. Time Post Treatment', 'FontSize', 14)
hold off

%% Displays Fit for 3 example weeks on one graph
v_model2_3 = model2pop16p( X, beta2newall, Vmaxbyweek, 3);
v_model2_6 = model2pop16p( X, beta2newall, Vmaxbyweek, 6);
v_model2_10 = model2pop16p( X, beta2newall, Vmaxbyweek, 10);
time(1, :) = 1:12;

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
%title ('Two Population Model 3, 6 & 10 Weeks Post Treatment','FontSize',18)
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
paramslb = zeros([16 1]);
paramsub = [ Inf; Inf; Inf; Inf; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1];

% sets up initial values of beta2new
params0 = [ .1; 17; .04; 50; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5];
%  res2sq=ones([n 1]);
[beta2new, resnorm2, residuals2] = lsqnonlin(@fit_simp2popabsresdensnormed,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    doseCV,...
    varCV,...
    k2,...
    wknumCV,...
    VmaxallCV);

beta2newCVtraining(k,:) = beta2new;
 % makes k by 16 parameter sets based on each training model
 [ lowerlimCV(:,k), upperlimCV(:,k), errorbarlengthCV(:,k), beta2newsimCV(:,k) ] = MCerrorinparams2( doseCV, beta2new, VmaxbyweekCV, k2, kstd, wknumCV, VmaxallCV);
 
v_model2trainingsetcell{k} =model2popallweeksnormed( doseCV, beta2new,VmaxbyweekCV, nsizeCV );
% uses each training set data and model paramters to create viability
% vector which represents that training models prediction
v_model2trainingsetplotcell{k} =model2popallweeksnormed( X, beta2new, VmaxbyweekCV,nsizeCV );
% does the same but across all doses given by X
end

coeffs2= beta2newCVtraining(1:num_cohorts, 1:4);
% k by 4 matrix of slope and center paramters from each training set model
% for clarity from fractional estimate parameters
%% Training Set Model Stand-Alone Statistics

% This section computes the validity of each training set model alone.  I
% ended up not including this in any analysis but still think it is
% indicative of the strengths of the different training models

for k = 1:num_cohorts
training_vddata = cell_of_trainingsets(k);
training_vddata = cell2mat(training_vddata);

doseCV = training_vddata(:,2);
varCV = training_vddata(:,3);

v_model2trainingset = v_model2trainingsetcell(k);
v_model2trainingset = cell2mat(v_model2trainingset);

actualresiduals2cell{k} = v_model2trainingset - varCV;
end

for k = 1:num_cohorts
actualresiduals2CV = actualresiduals2cell(k);
actualresiduals2CV = cell2mat(actualresiduals2CV);
nCV = length(actualresiduals2CV);
v_model2trainingset = v_model2trainingsetcell(k);
v_model2trainingset = cell2mat(v_model2trainingset);

sigma2CV(k) = std(actualresiduals2CV);
error2CV(k) = sum(abs(actualresiduals2CV))./nCV;
chi_squared2CV(k) = sum((actualresiduals2CV.^2)./v_model2trainingset);
RSS2(k) = sum(actualresiduals2CV.^2);
p2 = 16;
AIC2CV(k) = nCV*log(RSS2(k)./nCV) + 2*p2;
end


minAIC2CV = min(AIC2CV);
 % This next part was one way I had thought to develop model (by weighting
 % AIC values of each training set model) but was later told this was not
 % an accurate measure, so can be disregarded
 
for k = 1:num_cohorts
deltaAIC(k) = AIC2CV(k)- minAIC2CV;
end

for k = 1:num_cohorts
numweightsAIC2(k) = (exp(-0.5.*deltaAIC(k)));
end

for k = 1:num_cohorts
weightsAIC2(k) = numweightsAIC2(k)./(sum(numweightsAIC2));
model_contrib2(1,k) = weightsAIC2(k).*beta2newCVtraining(k,1);
model_contrib2(2,k) = weightsAIC2(k).*beta2newCVtraining(k,2);
model_contrib2(3,k) = weightsAIC2(k).*beta2newCVtraining(k,3);
model_contrib2(4,k) = weightsAIC2(k).*beta2newCVtraining(k,4);
model_contrib2(5,k) = weightsAIC2(k).*beta2newCVtraining(k,5);
model_contrib2(6,k) = weightsAIC2(k).*beta2newCVtraining(k,6);
model_contrib2(7,k) = weightsAIC2(k).*beta2newCVtraining(k,7);
model_contrib2(8,k) = weightsAIC2(k).*beta2newCVtraining(k,8);
model_contrib2(9,k) = weightsAIC2(k).*beta2newCVtraining(k,9);
model_contrib2(10,k) = weightsAIC2(k).*beta2newCVtraining(k,10);
model_contrib2(11,k) = weightsAIC2(k).*beta2newCVtraining(k,11);
model_contrib2(12,k) = weightsAIC2(k).*beta2newCVtraining(k,12);
model_contrib2(13,k) = weightsAIC2(k).*beta2newCVtraining(k,13);
model_contrib2(14,k) = weightsAIC2(k).*beta2newCVtraining(k,14);
model_contrib2(15,k) = weightsAIC2(k).*beta2newCVtraining(k,15);
model_contrib2(16,k) = weightsAIC2(k).*beta2newCVtraining(k,16);
end

beta2newCVfinal(1) = sum(model_contrib2(1,:));
beta2newCVfinal(2) = sum(model_contrib2(2,:));
beta2newCVfinal(3) = sum(model_contrib2(3,:));
beta2newCVfinal(4) = sum(model_contrib2(4,:));

fraction_estimates_sens(1) = sum(model_contrib2(5,:));
fraction_estimates_sens(2) = sum(model_contrib2(6,:));
fraction_estimates_sens(3) = sum(model_contrib2(7,:));
fraction_estimates_sens(4) = sum(model_contrib2(8,:));
fraction_estimates_sens(5) = sum(model_contrib2(9,:));
fraction_estimates_sens(6) = sum(model_contrib2(10,:));
fraction_estimates_sens(7) = sum(model_contrib2(11,:));
fraction_estimates_sens(8) = sum(model_contrib2(12,:));
fraction_estimates_sens(9) = sum(model_contrib2(13,:));
fraction_estimates_sens(10) = sum(model_contrib2(14,:));
fraction_estimates_sens(11) = sum(model_contrib2(15,:));
fraction_estimates_sens(12) = sum(model_contrib2(16,:));

params2newAIC = horzcat(beta2newCVfinal,fraction_estimates_sens);

fraction_estimates_res = ones([1 12]);
fraction_estimates_res = fraction_estimates_res - fraction_estimates_sens;
time = 1:1:12;

v_model2finalplot =model2popallweeksnormed( X, params2newAIC, Vmaxbyweek,nsize );
%% Calculate Average Training Sets

% Uses traing model paramters to come up with 95 percent confidence
% intervals of fitted parameters

beta2avgs = mean(beta2newCVtraining);
sigma2beta2 = std(beta2newCVtraining);
beta2avgs = beta2avgs';
sigma2beta2 = sigma2beta2';
beta295CI(:,1) = beta2avgs - 1.96.*sigma2beta2;
beta295CI(:,2) = beta2avgs + 1.96.*sigma2beta2;

sensfractsavg = beta2avgs(5:16);
resfractsavg = 1- beta2avgs(5:16);
resfracts95CI(:,1) = resfractsavg - 1.96.*sigma2beta2(5:16);
resfracts95CI(:,2) = resfractsavg + 1.96.*sigma2beta2(5:16);

figure(4)
hold off
plot(1:12, beta295CI(5:16, 1),'-r', 'LineWidth', 4)
hold on
plot(1:12, beta295CI(5:16, 2),'-r', 'LineWidth', 4)
plot (1:12, sensfractsavg, '-g', 'LineWidth', 4)
xlim ([1 12])
ylim ([ 0 1])
xlabel ('Weeks Post Treatment', 'FontSize', 20)
ylabel ('Fraction of Cells Sensitive', 'FontSize', 20)
%title ('95 Percent Confidence Interval of Sensitive Fraction','FontSize',20)

figure(5)
hold off
plot(1:12, resfracts95CI(:,1),'-r', 'LineWidth', 4)
hold on
plot(1:12, resfracts95CI(:,2),'-r', 'LineWidth', 4)
plot (1:12, resfractsavg, '-b', 'LineWidth', 4)
xlim ([1 12])
ylim ([ 0 1])
xlabel ('Weeks Post Treatment', 'FontSize', 20)
ylabel ('Fraction of Cells Resistant', 'FontSize', 20)
%title ('95 Percent Confidence Interval of Resistant Fraction','FontSize',20)

figure(6)
hold off
plot(time, resfractsavg, 'b', 'LineWidth', 4)
%errorbar(time, fres2, errorbarlength2(5:16)./2, 'b', 'LineWidth', 3)
hold on
plot(time, sensfractsavg, 'g', 'LineWidth',4)
%errorbar(time, fsens2, errorbarlength2(5:16)./2, 'g', 'LineWidth',3)
%legend('resistant (LD50 = 22.4 +/- 1.95 uM)', 'sensitive (LD50 = 77.4.1 +/- 3.65 uM)')
xlim([ 1 12])
ylim([0 1])
xlabel('Weeks Post Treatment', 'FontSize', 24)
ylabel('Fraction of Cells', 'FontSize',24)
%title ('Fraction of Resistant Cells and Sensitive Cells vs. Time Post Treatment', 'FontSize', 14)
hold off

%% Two Population General Model Testing Sets
% Now compare to the fit of a single cohort given by the testing set
% For each k, extract testing set 
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
nCVt = length(dose); % finds number of data points total in testing set

options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
%define bounds for all parameters
            % m1 cen1 m2 cen 2 wk1sens wk1res etc...
paramslb = zeros([12 1]);
paramsub = [ 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1];

% sets up initial values of beta2new
params0 = [ .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5];
%  res2sq=ones([n 1]);

[beta2new12p, resnorm2, residuals2] = lsqnonlin(@fit_simp2pop12pnormfix,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    doseCVt,...
    varCVt,...
    k2,...
    wknumCVt,...
    VmaxallCVt,...
    beta2avgs(1:4));
%coeffs2(k,:))
beta2newCVtesting(k,:) = beta2new12p;

end

removedgapsb2new = beta2newCVtesting;
timeind = removedgapsb2new == 0.5;
removedgapsb2new(timeind) = 0;
time = 1:1:12;
time = time';
for k = 1:num_cohorts
    fractionvt(:,1,k) = time;
    fractionvt(:,2,k)= removedgapsb2new(k,:)';
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

%% Take first training set and testing set and compute difference in model prediction

% Note, do not include in grant proposal code
% Must find the MC estimate of error in the parameters for each training
% set model. Then, can examine each testing set compared to the training
% set model to see if it falls within the CI

for k = 1:num_cohorts
beta2newtraining = beta2newCVtraining(k,:);

% This gives you a 3 column matrix (week, dose, viability) containing the model fit for every dose
% and every time point, now need to just find matching data points from
% corresponding test set

testingset = cell_of_testingsets(k);
testingset = cell2mat(testingset);

[ Vmaxbyweekt, Vmaxweekavgt, ninweekt, wknumt, Vmaxallt] = findVmaxandsizetest(testingset);

% puts training set parameters into testing set data
trsetpred = model2pop16ptesting( testingset, beta2newtraining, Vmaxallt );

testingset = horzcat(testingset, trsetpred(:,3)); % now testing set is 5 columns consisting of
 % wk, dose, expt viab, cohort, model fit
 
 residualstrvt = testingset(:,5) - testingset(:,3);
 testingset = horzcat(testingset,residualstrvt); 
 cell_to_compare{k} = testingset;
 
 chi_squaredtrvt(k) = sum((residualstrvt.^2)./testingset(:,5));
 
 figure(2*k-1)
% plot3(testingset(:,1), testingset(:,5), testingset(:,2), 'bo', 'LineWidth',4)
% hold on
% plot3(testingset(:,1), testingset(:,2), testingset(:,3), 'go')
 plot (testingset(:,2), testingset(:,5),'bo', 'LineWidth', 4)
 hold on
 plot (testingset(:,2), testingset(:,3), 'go')
% xlabel('time','FontSize',20) 
xlabel('dose (uM)', 'FontSize', 20)
ylabel('Viability', 'FontSize',20)
% title('Testing Set vs. Training Model', 'FontSize',20)
 
 figure(2*k)
 hold off
 plot(testingset(:,2), testingset(:,6), 'ro', 'LineWidth',4)
 ylim([-.5,.5])
end

% Now look at individual cohorts fit compared to training set model

figure(13)
plot3(testingset(:,1), testingset(:,2), testingset(:,5), 'bo', 'LineWidth',4)
hold on
plot3(testingset(:,1), testingset(:,2), testingset(:,3), 'go')

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
%title('Cohort Variation in Resistant Fraction', 'FontSize', 18)

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
%title('Cohort Variation in Sensitive Fraction', 'FontSize', 18)


%This returns a vector with each cohort in the z direction and each column
%is wpt, fsens, fres and each row is the value corresponding to that column
%% Determine initial subpopulation fractions


vddata_naive = load('naive_cohorts8910.m');

naive_cohorts = unique(vddata_naive(:,3));% gives cohort numbers (3,4,6,8,9)
num_naive_cohorts = length(naive_cohorts); % tells number of cohorts


%%
% for i = 1: num_naive_cohorts
%     
%     data_maniptest = vddata_naive;
%     index1 = data_maniptest(:,3) == naive_cohorts(i);
%     cohort_naive = data_maniptest(index1,:);
%     dosen = cohort_naive(:,1);
%     varn = cohort_naive(:,2);

   dosen8 = vddata_naive(1:12,1);
    varn8 = vddata_naive(1:12,2);
    
    dosen9 = vddata_naive(2:24,1);
    varn9 = vddata_naive(2:24,2);
    
    dosen10 = vddata_naive(25:36,1);
    varn10 = vddata_naive(25:36,2);

    
options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
%define bounds for all parameters
            % m1 cen1 m2 cen 2 wk1sens wk1res etc...
paramsnlb = 0;
paramsnub = 1;

% sets up initial values of beta2new
params0n = .5;


[betanaive, resnorm2, residuals2] = lsqnonlin(@fit_simp2popnaive,...
    params0n,...
    paramsnlb,...
    paramsnub,...
    options,...
    dosen10,...
    varn10,...
    beta2newCVfinal);

%betanaiveall(k,1) = betanaive;
