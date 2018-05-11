% This script makes Figure 1 (multiple options) for manuscript- an example of an exceptional
% cohort containing 12 data points at all time points for up to 15 weeks.

% Handles with graphics, 

clear all, close all, clc


vddata = load ('all_weeks_data_untreated.m'); % load concatenated data set (long Mx2 matrix
% all weeks added vertically, first column is week number
naive_vddata = load('naive_data.m');

%load cohort_8

time(1, :) = 1:15;
cohort_8ind = vddata(:,4) == 8;
cohort_8 = vddata(cohort_8ind, :);
% cohort_8 = load('cohort_8.mat');
% cohort_8 = struct2cell(cohort_8);
% cohort_8 = cell2mat(cohort_8); % Upload 203 x 4 matrix first column wk, 
%                                % second column dose, 3rd column viability,
%                                % 4th colum
cohort_8_naive = load('cohort_8_naive.mat');
cohort_8_naive = struct2cell(cohort_8_naive);
cohort_8_naive = cell2mat(cohort_8_naive);                          
 %% Allocate data
 dose = cohort_8(:,2);
var = cohort_8(:,3);
cohort_number = cohort_8(:,4);
n = length(dose); % finds number of data points total
wk = cohort_8(:,1);
[ Vmaxbyweek, Vmaxweekavg, ninweek, wknum, Vmaxall] = findVmaxandsizetest15(cohort_8);
nsize = wknum(:,2);
[ Vmaxbyweeknaive, Vmaxweekavgnaive, ninweeknaive, wknumnaive, Vmaxallnaive] = findVmaxandsizetest15(cohort_8_naive);
[kstd, Wstd, std_dev_model_fun] = fittechnicalerrors( dose);

dosenc = cohort_8_naive(:,2);
varnc = cohort_8_naive(:,3);

dose0ind = cohort_8_naive(:,2) == 0;
 Vmaxcohortnaivedata = cohort_8_naive(dose0ind,:);
 Vmaxnaiveavgcohort = mean(Vmaxcohortnaivedata(:,3));
 %% Find LD50 Naive

  options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
paramslbn = zeros([2 1]);
paramsubn = Inf([2 1]);

params0n = [ .1; 30];
% outputs betaLD50 with first 12 paramters slope values, second 12 are
% center values in week order
[betaLD50naive_8, resnormLD50naive_8, residualsLD50snaive_8] = lsqnonlin(@fit_simp1popnaive,...
    params0n,...
    paramslbn,...
    paramsubn,...
    options,...
    dosenc,...
    varnc,...
    Vmaxnaiveavgcohort);
 %% Fit cohort 8 to Time Constant Single Population Model
 
 initials1new = [ .05; 30]; % sets up initial values of beta1new

% call to lsqnonlin, which calls fit_simp1popnew which contains a
% difference vector in which the Vmax is replaced by different variables
% corresponding to the week. Here this parameter should be 1 every time
% (but its not)

[beta1_8unt, resnorm1_8, residuals1_8, output1, lambda1, Jacobian1] = lsqnonlin(@fit_simp1popunw, initials1new,[0; 0],[ Inf; Inf],[],dose,var, wknum, Vmaxall);
X = (0:.1:250);
X = X';
v_model1=model1pop( dose, beta1_8unt, Vmaxweekavg);
v_model1_naive = model1pop(dosenc, betaLD50naive_8, Vmaxnaiveavgcohort);
v_model1plot = model1pop(X, beta1_8unt, Vmaxweekavg);
actualresiduals1 = var-v_model1;
sigma1 = std(actualresiduals1);
error1 = sum(abs(actualresiduals1))./n;
chi_squared1 = sum((actualresiduals1.^2)./v_model1)
RSS1 = sum(actualresiduals1.^2);
p1 = 2;
AIC1 = n*log(RSS1./n) + 2*p1;
BIC1 = -n*log(RSS1./n) + p1*log(n);
DOF1 = n-p1;
MSE1 = chi_squared1./ DOF1;


%% Perform Monte Carlo Estimate of Error in Parameters
% based on standard deviation in technical replicates
[ lowerlim1MC, upperlim1MC, errorbarlength1MC, beta1simMC ] = MCerrorinparams1( cohort_8, v_model1, beta1_8, kstd, wknum, Vmaxall);
[lowerlimMCnaive, upperlimMCnaive, errorbarlengthMCnaive, betanaivesim] = MCerrorinparamsnaive( cohort_8_naive, v_model1_naive, betaLD50naive_8, kstd);
%% Print MC limits
lowerlim1MC = lowerlim1MC
upperlim1MC = upperlim1MC
lowerlimMCnaive = lowerlimMCnaive
upperlimMCnaive = upperlimMCnaive
%% Perform Boot Strapping Estimates of Error in Parameters
[lowerlim1boot, upperlim1boot] = BSerrorinparams1(residuals1_8, dose, beta1_8unt, Vmaxweekavg, cohort_8);
[lowerlim1bootnaive, upperlim1bootnaive] = BSerrorinparams1naive(residualsLD50snaive_8, dosenc, betaLD50naive_8, cohort_8_naive, Vmaxnaiveavgcohort);
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
title ('Single Population Model Cohort 8','FontSize',18)

figure(2)
hold off
plot(dose, residuals1_8,'o')
xlabel('dose (uM)','FontSize',18)
ylabel('Residuals','FontSize',18)
title ('Residuals for Single Population Model Cohort 8','FontSize',18)

 %% Fit cohort 8 to LD50s (Dynamic Single Population Model)
 
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
[betaLD50_8unt, resnormLD50_8, residualsLD50s_8] = lsqnonlin(@fit_simpLD5015,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    dose,...
    var,...
    wknum,...
    Vmaxall);

bulk_model_LD50_8 = model_bulk_LD50( dose, betaLD50_8unt, Vmaxall, wk );
%% Compute Statistics LD50(t)
actualresidualsLD50_8 = bulk_model_LD50_8-var;
sigmaLD50 = std(residualsLD50s_8)
errorLD50 = sum(abs(residualsLD50s_8))./n
chi_squaredLD50 = sum((residualsLD50s_8.^2)./bulk_model_LD50_8)
RSSLD50 = sum(residualsLD50s_8.^2);
pLD50 = 30;
AICLD50 = n*log(RSSLD50./n) + 2*pLD50
BICLD50 = -n*log(RSSLD50./n) + pLD50*log(n)
DOFLD50 = n-pLD50;
MSELD50 = chi_squaredLD50./ DOFLD50

%% Perform Monte Carlo Estimate of Error in Parameters LD50 (t)
% based on standard deviation in technical replicates
[ lowerlimLD50MC, upperlimLD50MC, errorbarlengthLD50MC, betaLD50simMC ] = MCerrorinparamsLD50( cohort_8, bulk_model_LD50_8, betaLD50_8, kstd, wknum, Vmaxall);

%% Print MC limits
lowerlimLD50MC = lowerlimLD50MC
upperlimLD50MC = upperlimLD50MC
errorbarlengthLD50MC = errorbarlengthLD50MC
%% Perform Boot Strapping Estimates of Error in Parameters
[lowerlimLD50boot, upperlimLD50boot] = BSerrorinparamsLD50(residualsLD50s_8, cohort_8, betaLD50_8unt);

%% Print Boot strapping limits
lowerlimLD50boot8 = lowerlimLD50boot
upperlimLD50boot8 = upperlimLD50boot
errorbarlengthLD50boot8untreated = upperlimLD50boot-lowerlimLD50boot

%% Plot of LD50s over time with bootstrapping error
figure (3)
errorbar(time, betaLD50_8unt(16:30),errorbarlengthLD50boot8untreated(16:30)./2, 'bo-', 'LineWidth',2)
hold on
errorbar(0, betaLD50naive_8(2), errorbarlength1bootnaive(2),'bo', 'LineWidth',2)
xlabel ('Time Post Treatment', 'FontSize', 20)
ylabel ('Single Population LD50', 'FontSize', 20)
ylim ([20,100])
title('Single Dynamic Population Model Cohort 8','FontSize', 20)

% Need to vertcat X 15 times
bulk_model_LD50_8_plot = model_bulk_LD50_plot( X, betaLD50_8, Vmaxbyweek);
X = (0:.1:250);
X = X';
X3D = vertcat(X,X,X,X,X,X,X,X,X,X,X,X,X,X,X);
t = X./(250./15);
t3D = vertcat(t,t,t,t,t,t,t,t,t,t,t,t,t,t,t);

figure(5)
hold off
plot3(X3D,t3D, bulk_model_LD50_8_plot, 'o'); 
xlabel(' Dose Doxorubicin (uM)', 'FontSize', 12)
zlabel('Model Prediction and Experiment Viability', 'FontSize',12)
ylabel('Time (Weeks Post Treatment)','FontSize', 12)

figure(4)
hold off
plot3(dose, cohort_8(:,1), bulk_model_LD50_8, 'ko', 'LineWidth',5)
hold on
plot3(dose, cohort_8(:,1), var, 'bo')
xlabel(' Dose Doxorubicin (uM)', 'FontSize', 12)
zlabel('Model Prediction and Experiment Viability', 'FontSize',12)
ylabel('Time (Weeks Post Treatment)','FontSize', 12)

%% Fit cohort 8 to two population model

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
[beta2_8unt, resnorm2_8, residuals2_8] = lsqnonlin(@fit_simp2popunw15,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    dose,...
    var,...
    wknum,...
    Vmaxall);

%% Compute and Output Two Population Statistics and Model Comparison

v_model2allweeksnormed_8 = model2popallweeksnormed( dose, beta2_8unt, Vmaxbyweek, nsize);
actualresiduals2 = var- v_model2allweeksnormed_8;
sigma2 = std(residuals2_8)
error2 = sum(abs(residuals2_8))./n
chi_squared2 = sum((residuals2_8.^2)./v_model2allweeksnormed_8)
RSS2 = sum(residuals2_8.^2);
p2 = 19;
AIC2 = n*log(RSS2./n) + 2*p2
BIC2 = -n*log(RSS2./n) + p2*log(n)

DOF2 = n-19;
MSE2 = chi_squared2./(DOF2)
%% Find Naive Overall Fractions

% simply want to average all data to find initial starting conditions of
% all cohorts
params0 = .5;
paramslb = 0;
paramsub = 1;
coeffs2 = beta2_8(1:4);
[beta2naive_8, resnorm2naive_8, residuals2naive_8] = lsqnonlin(@fit_simp2popnaiveunw,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    dosenc,...
    varnc,...
    coeffs2);

fsens0 = beta2naive_8;
fres0 = 1-beta2naive_8;

v_model2naive = model2popallweeksnormednaive( dosenc, coeffs2, beta2naive_8, Vmaxnaiveavgcohort)
%%
fres2 = zeros([15,1]);
fsens2 = beta2_8(5:19);

for i = 1:15
    fres2(i) = 1-fsens2(i);
end
time = zeros([1 15]);
time(1, :) = 1:15;


%% Perform Monte Carlo Estimate of Error in Parameters two population model
% based on standard deviation in technical replicates
 [ lowerlim2MC, upperlim2MC, errorbarlength2MC, beta2simMC ] = MCerrorinparams2( cohort_8, v_model2allweeksnormed_8, beta2_8, Vmaxbyweek, kstd, wknum, Vmaxall);
 [ lowerlim2MCnaive, upperlim2MCnaive, errorbarlength2MCnaive, beta2naivesimMC ] = MCerrorinparams2naive( cohort_8_naive, v_model2naive, beta2_8, beta2naive_8, kstd);
%% Print MC limits two population
lowerlim2MC = lowerlim2MC
upperlim2MC = upperlim2MC
errorbarlength2MC = errorbarlength2MC
lowerlim2MCnaive = lowerlim2MCnaive
upperlim2MCnaive = upperlim2MCnaive
errorbarlength2MCnaive = errorbarlength2MCnaive

%% Perform Boot Strapping Estimates of Error in Parameters two population model
% Need to make this function and its naive counterpart
[lowerlim2boot8, upperlim2boot8] = BSerrorinparams2(residuals2_8, beta2_8unt, cohort_8);
[lowerlim2bootnaive8, upperlim2bootnaive8] = BSerrorinparams2naive(residuals2naive_8, dosenc, beta2naive_8, cohort_8_naive, beta2_8, Vmaxnaiveavgcohort);
%% Print Boot strapping limits
lowerlim2boot8 = lowerlim2boot8
upperlim2boot8 = upperlim2boot8
errorbarlength2boot8unt = upperlim2boot8-lowerlim2boot8
lowerlim2bootnaive8 = lowerlim2bootnaive8
upperlim2bootnaive8 = upperlim2bootnaive8
errorbarlength2bootnaive8 = upperlim2bootnaive8 - lowerlim2bootnaive8
%% Plots of sensitive and resistant fractions over time with bootstrapping error



% Plot of the resistant and sensitive fractions
fig6 = figure;
fig7 = figure;
figure(fig6);
hold off
errorbar(time, fres2, errorbarlength2boot8(5:19)./2, 'bo-', 'LineWidth', 2);
hold on
errorbar(0, fres0, errorbarlength2bootnaive8./2,'bo','LineWidth',5)
%errorbar(time, fsens2, errorbarlength2boot8(5:19)./2, 'go-', 'LineWidth', 2)
%errorbar(0, fsens0, errorbarlength2bootnaive8./2,'go','LineWidth',5)
xlim([ 0 15])
ylim([ 0 1])
xlabel('Time (Weeks Post Treatment)', 'FontSize', 24)
ylabel('Fraction of Cells', 'FontSize',24)
%title ('Fraction of Resistant Cells and Sensitive Cells vs. Time Post Treatment from All Data', 'FontSize', 14)
hold off


figure(fig7);
hold off
errorbar(time, fsens2, errorbarlength2boot8(5:19)./2, 'go-', 'LineWidth', 2)
hold on
errorbar(0, fsens0, errorbarlength2bootnaive8./2,'go','LineWidth',5)
xlim([ 0 15])
ylim([ 0 1])
xlabel('Time (Weeks Post Treatment)', 'FontSize', 24)
ylabel('Fraction of Cells', 'FontSize',24)
%title ('Fraction of Resistant Cells and Sensitive Cells vs. Time Post Treatment from All Data', 'FontSize', 14)
hold off
%%
get(fig6)
set(fig6)
c = [1.0 0.8 0.0]
set(fig6, 'Color','w'); % changes background color
% can change the name of figure
% can set position via
fig6 = figure('Color','w', 'Position', [494 87 560 420])
%subplot(rows, columns, number)
%%
fig1 = figure('Color','w');
h = errorbar(time, fsens2, errorbarlength2boot8(5:19)./2, 'go-', 'LineWidth', 2)
hold on
h = errorbar(0, fsens0, errorbarlength2bootnaive8./2,'go','LineWidth',2)
h2 = get (h, 'Parent')
set(h2, 'FontSize',14,'LineWidth',2);
xlabel('Time (Weeks Post Treatment)')
ylabel('Fraction of Sensitive Cells')
title ('Two Population Model Cohort 8: Sensitive Fraction')
xlim([0 16])
ylim([0 1])
%%
fig2 = figure('Color','w');
h = errorbar(time, fres2, errorbarlength2boot8(5:19)./2, 'bo-', 'LineWidth', 2)
hold on
h = errorbar(0, fres0, errorbarlength2bootnaive8./2,'bo','LineWidth',2)
h2 = get (h, 'Parent')
set(h2, 'FontSize',14,'LineWidth',2);
xlabel('Time (Weeks Post Treatment)')
ylabel('Fraction of Resistant Cells')
title ('Two Population Model Cohort 8: Resistant Fraction')
xlim([0 16])
ylim([0 1])

