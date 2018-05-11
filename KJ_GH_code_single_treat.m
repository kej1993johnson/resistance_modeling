% This script loads in data that is copy and pasted from excel into a
% matlab file. The data is in the format of:
% first column: time post treatment (weeks): wk
% second column: dose administered (nM): dose
% third column: fraction viable at that dose: var
% fourth column: cohort (look back at excel to see dates): cohort_number

% This script loads either the data containing cohorts 3,4, 6,8,9, and 10
% or cohorts 9, 10, and 11 in response to a single treatment at 0 weeks.
% The analysis goes out to 8 weeks.


clear all, close all, clc
%vddata = load ('eight_weeks_data_2_18.m'); % Data from cohorts 3,4,6,8,9,10 
%(used to make figures in manuscript draft/on poster
vddata = load ('single_treat_data_cohorts91011.m'); % data from cohorts 9, 10 and 11 
% use this data to make figures of single treatment
% all weeks added vertically, first column is week number
naive_vddata = load('naive_data.m'); % loads the 0th week dose-response for cohorts 3,4,6,8,9,10
%% Remove Cohort 3
% Run this if want to remove the third cohort from the analysis of cohorts
% 3,6,8,9 and 10. Otherwise do not
coh_ind = vddata(:, 4) == 3;
vddata = vddata(coh_ind == 0, :);

%% Allocate Variables
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


%% Find LD50 of Week 0 (Naive) Cells

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

[beta1new, resnorm1, residuals1] = lsqnonlin(@fit_simp1popunw, initials1new,[0; 0],[ Inf; Inf],[],dose,var, wknum, Vmaxall);

% Plugs parameters back into model to plot then calculates model statistics
X = (0:.1:250);
X = X';
v_model1=model1pop( dose, beta1new, Vmaxweekavg );
v_model1plot = model1pop(X, beta1new, Vmaxweekavg);
actualresiduals1 = v_model1-var;
sigma1 = std(residuals1);
error1 = sum(abs(residuals1))./n;
chi_squared1 = sum((residuals1.^2)./v_model1);
RSS1 = sum(residuals1.^2);
p1 = 2;
AIC1 = n*log(RSS1./n) + 2*p1;
BIC1 = -n*log(RSS1./n) + p1*log(n);
DOF1 = n-p1;
MSE1 = chi_squared1./ DOF1;

plot(dose, var, '.', X, v_model1plot);

%% Perform Boot Strapping Estimates of Error in Parameters for Single LD50 model and for week 0 LD50
[lowerlim1boot, upperlim1boot] = BSerrorinparams1(residuals1, dose, beta1new, Vmaxweekavg, vddata);
[lowerlim1bootnaive, upperlim1bootnaive] = BSerrorinparams1naive(residualsLD50snaive, dosenc, betaLD50naive, naive_vddata, Vmaxnaiveavg);
errorbarlength1boot = upperlim1boot-lowerlim1boot;
errorbarlength1bootnaive = upperlim1bootnaive-lowerlim1bootnaive;


%% Plot Bulk LD50

% Normally want to compare this to the untreated cohorts at each time, so
% we load those LD50s fro the script ("fit2modelscript8weeksuntreated")
% where these are saved as beta and var

% This will produce an error when trying for only cohorts 9,10, and 11
% since the naive data is made up of all of the cohor
beta1unt = load('beta1unt_8wks.mat');
beta1unt = struct2cell(beta1unt);
beta1unt = cell2mat(beta1unt);

dose_unt = load('dose_unt_8wks.mat');
dose_unt = struct2cell(dose_unt);
dose_unt = cell2mat(dose_unt);

var_unt = load('var_unt_8wks.mat');
var_unt = struct2cell(var_unt);
var_unt = cell2mat(var_unt);

v_model1plotunt = model1pop(X, beta1unt, Vmaxweekavg);
figure(1) % graphs fit to all weeks in red
hold off
plot (dose, var, 'b.', 'LineWidth', 1.5)
hold on
plot(dose_unt,var_unt, 'k.', 'LineWidth', 1.5)
plot (X, v_model1plot, '-b', 'LineWidth', 3);
plot(X, v_model1plotunt, '-k', 'LineWidth',3);
xlim([0,250]);
xlabel('Dose (uM)','FontSize',18)
ylabel('Viability','FontSize',18)
title ('Single Static Population Model','FontSize',18)
legend('Experimental Treated Viability','Experimental Untreated Viability','Treated Single Static Population Model', 'Untreated Single Static Population Model')
set(gca,'LineWidth',1.5,'FontSize',12)


%% Dynamic Single Population Model: 8 weeks of data

options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
%define bounds for all parameters
            % m1 cen1 m2 cen 2 wk1sens wk1res etc...
paramslb = zeros([16 1]);
paramsub = Inf([16 1]);


params0 = [ .1; .1; .1; .1; .1; .1; .1; .1; 30; 30; 30; 30; 30; 30; 30; 30];
% outputs betaLD50 with first 12 paramters slope values, second 12 are
% center values in week order
[betaLD50, resnormLD50, residualsLD50] = lsqnonlin(@fit_simpLD508,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    dose,...
    var,...
    wknum,...
    Vmaxall);

%% Perform Boot Strapping Estimates of Error in Parameters Single Dynamic Model
[lowerlimLD50boot, upperlimLD50boot] = BSerrorinparamsLD50(residualsLD50, vddata, betaLD50);

lowerlimLD50boot = lowerlimLD50boot;
upperlimLD50boot = upperlimLD50boot;
errorbarlengthLD50boot = upperlimLD50boot-lowerlimLD50boot;



%% Plot of LD50s over time
betaLD50unt = load('betaLD50unt_8wks.mat');
betaLD50unt = struct2cell(betaLD50unt);
betaLD50unt = cell2mat(betaLD50unt);

errorbarlengthLD50untreated = load('errorbarlengthLD50bootuntreated_8wks.mat');
errorbarlengthLD50untreated = struct2cell(errorbarlengthLD50untreated);
errorbarlengthLD50untreated = cell2mat(errorbarlengthLD50untreated);

errorbarsLD50 = horzcat(errorbarlength1bootnaive(2),errorbarlengthLD50boot(9:16));
errorbarsLD50unt = horzcat(errorbarlength1bootnaive(2), errorbarlengthLD50untreated(9:16));
time = 0:1:8;
LD50s = vertcat(betaLD50naive(2), betaLD50(9:16));
LD50sunt = vertcat(betaLD50naive(2), betaLD50unt(9:16));

figure (2)
set(gca,'LineWidth',1.5,'FontSize',12);
hold off
% Plot of LD50s over time
%plot(time, betaLD50(16:30), 'o', 'LineWidth', 4)
errorbar(time, LD50s, errorbarsLD50./2, '-go', 'LineWidth', 2)
hold on
%errorbar(0, betaLD50naive(2), errorbarlength1bootnaive(2)./2,'ko', 'LineWidth',2)
%plot(0, betaLD50naive(2), 'go', 'LineWidth',4)
xlabel ('Time post treatment (weeks)', 'FontSize', 16)
ylabel ('Single Population LD50 (uM)', 'FontSize', 16)
set(gca,'LineWidth',1.5,'FontSize',16)
errorbar(time, LD50sunt, errorbarsLD50unt./2,'ko', 'LineWidth',2)
ylim ([20,90])
title('Single Dynamic Population Model','FontSize', 18)
legend('Treated Singe Population LD50', 'Untreated Single Population LD50')




%% Model Fit Single Population Shifting in Time
bulk_model_LD50 = model1popallweeksnormed( dose, betaLD50, Vmaxbyweek, nsize);
% makes curve for each week

% Plot example model fits of single dynamic model
time = 1:1:8;


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
% Model statistics single dynamic to single static
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

plot(dose, var, '.', dose, bulk_model_LD50, 'r.')

% Comparison of one population changing in time with constant 1 population model
F_stat1vLD50alt = ((RSS1-RSSLD50)./(DOF1-DOFLD50))./(RSSLD50./(DOFLD50))
p1vLD50alt = 1-fcdf(F_stat1vLD50alt, DOF1-DOFLD50, DOFLD50)


%% Two Population Model fractions shift in time

options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
%define bounds for all parameters
            % m1 cen1 m2 cen 2 wk1sens wk1res etc...
paramslb = zeros([12 1]);
paramsub = [ Inf; Inf; Inf; Inf; 1; 1; 1; 1; 1; 1; 1; 1];

% sets up initial values of beta2new
params0 = [ .1; 17; .04; 50; .5; .5; .5; .5; .5; .5; .5; .5];
%  res2sq=ones([n 1]);
[beta2newall, resnorm2, residuals2] = lsqnonlin(@fit_simp2popunw8,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    dose,...
    var,...
    wknum,...
    Vmaxall);
coeffs2 = beta2newall(1:4);
%% Compute and Output Two Population Statistics and Model Comparison

v_model2allweeksnormed = model2popallweeksnormed( dose, beta2newall, Vmaxbyweek, nsize);
% again plots the model prediction at each dose for all weeks
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
% Comparison of two population model with bulk constant time model
F_stat1v2alt = ((RSS1-RSS2)./(DOF1-DOF2))./(RSS2./(DOF2))
p1v2alt = 1-fcdf(F_stat1v2alt, DOF1-DOF2, DOF2)

% Comparison of one population changing in time with two population model
F_stat2vLD50alt = ((RSS2-RSSLD50)./(DOF2-DOFLD50))./(RSSLD50./(DOFLD50))
p2vLD50alt = 1-fcdf(F_stat2vLD50alt, DOF2-DOFLD50, DOFLD50)


%2D view
figure(10)
hold off
plot(dose, var, 'bo')
hold on
plot(dose, v_model2allweeksnormed,'ro')

%3D view
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

%% Find Naive (week 0) Fres and Fsens
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
    coeffs2,...
    Vmaxnaiveavg);

fsens0 = beta2naive;
fres0 = 1-beta2naive;
%Find Naive model
v_model2allweeksnormednaive = model2popallweeksnormednaive( dose, coeffs2, beta2naive)

%% Perform Boot Strapping Estimates of Error in Parameters two population model
% Need to make this function and its naive counterpart
[lowerlim2boot, upperlim2boot] = BSerrorinparams2(residuals2, beta2newall, vddata);
[lowerlim2bootnaive, upperlim2bootnaive] = BSerrorinparams2naive(residuals2naive, dosenc, beta2naive, naive_vddata, beta2newall, Vmaxnaiveavg);
lowerlim2boot = lowerlim2boot;
upperlim2boot = upperlim2boot;
errorbarlength2boot = upperlim2boot-lowerlim2boot;
lowerlim2bootnaive = lowerlim2bootnaive;
upperlim2bootnaive = upperlim2bootnaive;
errorbarlength2bootnaive = upperlim2bootnaive - lowerlim2bootnaive;

%% Plot of resistant and sensitive fraction estimates over time

beta2unt = load('beta2unt_8wks.mat');
beta2unt = struct2cell(beta2unt);
beta2unt = cell2mat(beta2unt);

errorbarlength2bootuntreated = load('errorbarlength2bootuntreated_8wks.mat');
errorbarlength2bootuntreated = struct2cell(errorbarlength2bootuntreated);
errorbarlength2bootuntreated = cell2mat(errorbarlength2bootuntreated);
fres2 = zeros([8,1]);
fsens2 = beta2newall(5:12);

for i = 1:8
    fres2(i) = 1-fsens2(i);
end
time = zeros([1 9]);
time(1, :) = 0:1:8;

fres2unt = zeros([8,1]);
fsens2unt = beta2unt;

for i = 1:8
    fres2unt(i) = 1-fsens2unt(i);
end
fres2 = vertcat(fres0, fres2);
fsens2 = vertcat(fsens0, fsens2);

fres2unt = vertcat(fres0, fres2unt);
fsens2unt = vertcat(fsens0, fsens2unt);

errorbars = horzcat(errorbarlength2bootnaive, errorbarlength2boot(5:12));
errorbars_unt = horzcat(errorbarlength2bootnaive, errorbarlength2bootuntreated(5:12));

% Plot of the resistant and sensitive fractions
subplot(1,2,1)
hold off
errorbar(time, fres2, errorbars./2, 'bo-', 'LineWidth', 1.5)
hold on
%errorbar(time, fsens2, errorbars./2, 'go-', 'LineWidth', 1.5)
errorbar(time, fres2unt, errorbars_unt./2, 'ko', 'LineWidth', 1.5)
%errorbar(0, fres0, errorbarlength2bootnaive./2, 'bo', 'LineWidth', 1.5)
%errorbar(time, fres2unt, errorbars_unt./2, 'ko-', 'LineWidth', 2)
%errorbar(0, fsens0, errorbarlength2bootnaive./2, 'ro', 'LineWidth', 1.5)
ylim([ 0 1])
xlabel('Time Post Treatment(Weeks)', 'FontSize', 14)
ylabel('Fraction of Cells', 'FontSize',14)
title ('Fraction of Resistant Cells Two Population Model', 'FontSize', 14)
set(gca,'LineWidth',1.5,'FontSize',12)
legend ('Fraction of Treated Cells Resistant', 'Fraction of Untreated Cells Resistant')
hold off

subplot(1,2,2)
hold off
hold on
errorbar(time, fsens2, errorbars./2, 'yo-', 'LineWidth', 1.5)
%errorbar(time, fres2unt, errorbars_unt./2, 'ko', 'LineWidth', 1.5)
%errorbar(0, fres0, errorbarlength2bootnaive./2, 'bo', 'LineWidth', 1.5)
errorbar(time, fsens2unt, errorbars_unt./2, 'ko', 'LineWidth', 2)
%errorbar(0, fsens0, errorbarlength2bootnaive./2, 'ro', 'LineWidth', 1.5)
ylim([ 0 1])
xlabel('Time Post Treatment (Weeks)', 'FontSize', 14)
ylabel('Fraction of Cells', 'FontSize',14)
title ('Fraction of Sensitive Cells Two Population Model', 'FontSize', 14)
set(gca,'LineWidth',1.5,'FontSize',12)
legend ('Fraction of Treated Cells Sensitive', 'Fraction of Untreated Cells Sensitive')
hold off
%% Model Fit for 2 example weeks on one graph
v_model2_2 = model2pop12p( X, beta2newall, 2, Vmaxbyweek);
v_model2_8 = model2pop12p( X, beta2newall,  8, Vmaxbyweek);
%v_model2_10 = model2pop16p( X, beta2newall, Vmaxbyweek, 10);


figure(7) % looks at example weeks 3, 6, and 10 WPT against all experimental data
hold off
plot (dose(sum(nsize(1:1))+1:sum(nsize(1:2))), var(sum(nsize(1:1))+1:sum(nsize(1:2))), 'ro')
hold on
%plot (dose(sum(nsize(1:5))+1:sum(nsize(1:6))), var(sum(nsize(1:5))+1:sum(nsize(1:6))), 'go')
%plot (dose(sum(nsize(1:9))+1:sum(nsize(1:10))), var(sum(nsize(1:9))+1:sum(nsize(1:10))), 'bo')
plot (X, v_model2_2, '-r', 'LineWidth', 3);
plot (X, v_model2_8, '-g', 'LineWidth', 3);
%plot (X, v_model2_10, '-b', 'LineWidth', 3);
xlim ([0 250])
xlabel('dose (uM)','FontSize',24)
ylabel('viability','FontSize',24)
title ('Two Population Model 3, 6 & 10 Weeks Post Treatment','FontSize',18)