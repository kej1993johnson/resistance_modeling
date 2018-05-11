% This script is for the second treatment at 6 weeks and does analysis out
% to 15 weeks 

clear all, close all, clc

%vddata = load ('secondtreatment.m'); % only dosed once ones
vddata = load ('second_treatment_15wks.m'); % dosed twice
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

%[new_vddata] = fit_by_var_by_time(vddata);
figure(10)
hold off
plot(dose, new_vddata(:,5),'o');
%ylim([0 100])
xlabel('dose(uM)', 'FontSize', 20)
ylabel('Weight by Variance','FontSize', 20)
title('Plot of Proposed Weighting Scheme', 'FontSize', 20)

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

[beta1new, resnorm1, residuals1] = lsqnonlin(@fit_simp1popunw, initials1new,[0; 0],[ Inf; Inf],[],dose,var, wknum, Vmaxall);
X = (0:.1:250);
X = X'
v_model1=model1pop( dose, beta1new, Vmaxweekavg );
v_model1plot = model1pop(X, beta1new, Vmaxweekavg);
actualresiduals1 = v_model1-var;
sigma1 = std(actualresiduals1)
error1 = sum(abs(actualresiduals1))./n
chi_squared1 = sum((actualresiduals1.^2)./v_model1)
RSS1 = sum(actualresiduals1.^2);
p1 = 2;
AIC1 = n*log(RSS1./n) + 2*p1
BIC1 = -n*log(RSS1./n) + p1*log(n)
DOF1 = n-p1;
MSE1 = chi_squared1./ DOF1

%% Perform Boot Strapping Estimates of Error in Parameters
[lowerlim1boot, upperlim1boot] = BSerrorinparams1(residuals1, dose, beta1new, Vmaxweekavg, vddata);
[lowerlim1bootnaive, upperlim1bootnaive] = BSerrorinparams1naive(residualsLD50snaive, dosenc, betaLD50naive, naive_vddata, Vmaxnaiveavg);
lowerlim1boot = lowerlim1boot;
upperlim1boot = upperlim1boot;
lowerlim1bootnaive = lowerlim1bootnaive;
upperlim1bootnaive = upperlim1bootnaive;
errorbarlength1boot = upperlim1boot-lowerlim1boot;
errorbarlength1bootnaive = upperlim1bootnaive-lowerlim1bootnaive;

%% Plot Bulk LD50

beta1unt = load('beta1unt.mat');
beta1unt = struct2cell(beta1unt);
beta1unt = cell2mat(beta1unt);

beta1onedose = load('beta1onedose.mat');
beta1onedose = struct2cell(beta1onedose);
beta1onedose = cell2mat(beta1onedose);


v_model1plotunt = model1pop(X, beta1unt, Vmaxweekavg);
v_model1plotonedose = model1pop(X, beta1onedose, Vmaxweekavg);
figure(1) % graphs fit to all weeks in red
hold off
plot (dose, var, 'b.')
hold on
plot (X, v_model1plot, '-r', 'LineWidth', 3);
plot(X, v_model1plotonedose, '-g', 'LineWidth',3);
plot(X, v_model1plotunt, '-k', 'LineWidth',3);
xlim([0,250]);
xlabel('Dose (uM)','FontSize',18)
ylabel('Viability','FontSize',18)
title ('Single Static Population Model Effect of Pulse Treatment','FontSize',14)
legend('Second Treatment Data','Twice Treated LD50 = 57.47 uM', 'Single Treatment LD50 = 51.30 uM', 'Untreated LD50 = 38.18 uM')
set(gca,'LineWidth',1.5,'FontSize',12)


%% Dynamic Single Population Model: 15 weeks of data

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
[betaLD50twodose, resnormLD50, residualsLD50] = lsqnonlin(@fit_simpLD5015,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    dose,...
    var,...
    wknum,...
    Vmaxall);

vmodel1 = model1popallweeksnormed( dose, betaLD50twodose, Vmaxbyweek, nsize);

plot(dose, vmodel1, 'o')

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
[lowerlimLD50boot, upperlimLD50boot] = BSerrorinparamsLD50(residualsLD50, vddata, betaLD50twodose);
lowerlimLD50boot = lowerlimLD50boot;
upperlimLD50boot = upperlimLD50boot;
errorbarlengthLD50boot = upperlimLD50boot-lowerlimLD50boot;


%% Plot of LD50s over time
betaLD50unt = load('betaLD50unt.mat');
betaLD50unt = struct2cell(betaLD50unt);
betaLD50unt = cell2mat(betaLD50unt);

errorbarlengthLD50untreated = load('errorbarlengthLD50bootuntreated.mat');
errorbarlengthLD50untreated = struct2cell(errorbarlengthLD50untreated);
errorbarlengthLD50untreated = cell2mat(errorbarlengthLD50untreated);

betaLD50onedose = load('betaLD50onedose.mat');
betaLD50onedose = struct2cell(betaLD50onedose);
betaLD50onedose = cell2mat(betaLD50onedose);

errorbarlengthLD50onedose = load('errorbarlengthLD50onedose.mat');
errorbarlengthLD50onedose = struct2cell(errorbarlengthLD50onedose);
errorbarlengthLD50onedose = cell2mat(errorbarlengthLD50onedose);

errorbarsLD50= errorbarlengthLD50boot';

errorbarsLD50onedose = errorbarlengthLD50onedose';

time = 0:1:15;
betaLD50twodose = vertcat(betaLD50naive(2), betaLD50twodose(16:30));
errorbarLD50twodose = vertcat(errorbarlength1bootnaive(2), errorbarsLD50(16:30));
%%
betaLD50unt = vertcat(betaLD50naive(2), betaLD50unt(16:30));
errorbarsLD50untreated = vertcat( errorbarlength1bootnaive(2), errorbarlengthLD50untreated(16:30)');

betaLD50onedose = vertcat(betaLD50naive(2), betaLD50onedose(16:30));
errobarsLD50onedose1 = vertcat(errorbarlength1bootnaive(2), errorbarsLD50onedose(16:30));

%%
figure (2)
hold off
% Plot of LD50s over time
%plot(time, betaLD50(16:30), 'o', 'LineWidth', 4)
errorbar(time, betaLD50twodose, errorbarLD50twodose./2, '-ro', 'LineWidth', 2)
hold on
errorbar(time, betaLD50onedose, errobarsLD50onedose1./2, '-go', 'LineWidth',2)
xlabel ('Time Post Treatment (Weeks)', 'FontSize', 16)
ylabel ('Single Population LD50', 'FontSize', 16)
set(gca,'LineWidth',1.5,'FontSize',14)
errorbar(time, betaLD50unt, errorbarsLD50untreated./2,'ko', 'LineWidth',2)
ylim ([20,85])
title('Single Dynamic Population Model Effect of Pulse Treatments ','FontSize', 18)
legend('Twice Treated at 0 and 6 weeks','Single Treatment 0 weeks', 'Untreated Response')

%% Look at response after first and second treatment side by side
figure (3)
hold off
% Plot of LD50s over time
%plot(time, betaLD50(16:30), 'o', 'LineWidth', 4)
errorbar(time(1:10), betaLD50twodose(7:16), errorbarLD50twodose(7:16)./2, '-ro', 'LineWidth', 2)
hold on
errorbar(time(1:10), betaLD50onedose(1:10), errobarsLD50onedose1(1:10)./2, '-go', 'LineWidth',2)
xlabel ('Time Post Treatment (Weeks)', 'FontSize', 16)
ylabel ('Single Population LD50', 'FontSize', 16)
set(gca,'LineWidth',1.5,'FontSize',14)
errorbar(time(1:10), betaLD50unt(1:10), errorbarsLD50untreated(1:10)./2,'ko', 'LineWidth',2)
ylim ([20,85])
title('Single Dynamic Population Model Response to First and Second Treatment ','FontSize', 16)
legend('Second Treatment Response','First Treatment Response', 'Untreated Response')


%% Compute Statistics: Comparing Stagnant One Population to Dynamic One Population

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



%% Two Dose Two population model
%One Dose
beta2onedose = load('beta2onedose.mat');
beta2onedose = struct2cell(beta2onedose);
beta2onedose = cell2mat(beta2onedose);

errorbarlength2onedose = load('errorbarlength2onedose.mat');
errorbarlength2onedose = struct2cell(errorbarlength2onedose);
errorbarlength2onedose = cell2mat(errorbarlength2onedose);
%%
coeffs2 = load('coeffs2.mat'); % Coefficients (Parameters) taken from single treatment model ('KJ_GH_code_single_treat.m)
coeffs2 = struct2cell(coeffs2);
coeffs2 = cell2mat(coeffs2);
options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
%define bounds for all parameters
            % m1 cen1 m2 cen 2 wk1sens wk1res etc...
paramslb = zeros([15 1]);
paramsub = [ Inf; Inf; Inf; Inf; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1];

% sets up initial values of beta2new
params0 = [ .1; 17; .04; 50; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5];
%  res2sq=ones([n 1]);
[beta2twodose, resnorm2, residuals2] = lsqnonlin(@fit_simp2popunw15_untreated,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    dose,...
    var,...
    wknum,...
    Vmaxall,...
    coeffs2);
%% BS error second dose two population
beta2_comb = vertcat(coeffs2, beta2twodose);
[lowerlim2boot, upperlim2boot] = BSerrorinparams2(residuals2, beta2_comb, vddata);
errorbarlength2twodose = upperlim2boot-lowerlim2boot;
%% Find Naive Overall Fractions
% loads 3 column vector of dose, viability, and cohort number
 
% simply want to average all data to find initial starting conditions of
% all cohorts
params0 = .5;
paramslb = 0;
paramsub = 1;
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

%% BS error naive
[lowerlim2bootnaive, upperlim2bootnaive] = BSerrorinparams2naive(residuals2naive, dosenc, beta2naive, naive_vddata, beta2onedose, Vmaxnaiveavg);
errorbarlength2bootnaive = upperlim2bootnaive - lowerlim2bootnaive
%% Compile untreated, single dose, and repeat dose coefficients
%Untreated
beta2unt = load('beta2unt.mat');
beta2unt = struct2cell(beta2unt);
beta2unt = cell2mat(beta2unt);

fres2unt = zeros([15,1]);
fsens2unt = beta2unt;

for i = 1:15
    fres2unt(i) = 1-fsens2unt(i);
end
fres2untall = vertcat(fres0,fres2unt);

fres2onedose = zeros([15,1]);
fsens2onedose = beta2onedose;

for i = 1:15
    fres2onedose(i) = 1-fsens2onedose(i);
end
fres2onedoseall = vertcat(fres0, fres2onedose);

fres2twodose = zeros([15,1]);
fsens2twodose = beta2twodose;

for i = 1:15
    fres2twodose(i) = 1-fsens2twodose(i);
end

fres2twodoseall = vertcat(fres0, fres2twodose);

errorbarlength2bootuntreated = load('errorbarlength2bootuntreated.mat');
errorbarlength2bootuntreated = struct2cell(errorbarlength2bootuntreated);
errorbarlength2bootuntreated = cell2mat(errorbarlength2bootuntreated);
errorbar2untreated = vertcat( errorbarlength2bootnaive, errorbarlength2bootuntreated(5:19)');
%errorbar2untreated is 16 weeks of fraction estimates

errorbar2onedose = vertcat(errorbarlength2bootnaive, errorbarlength2onedose(5:19)');
errorbar2twodose = vertcat(errorbarlength2bootnaive, errorbarlength2onedose(5:19)');


time = 0:1:15;
figure (4)
hold off
% Plot of LD50s over time
%plot(time, betaLD50(16:30), 'o', 'LineWidth', 4)
errorbar(time, fres2twodoseall, errorbar2twodose./2, '-ro', 'LineWidth', 2)
hold on
errorbar(time,fres2onedoseall, errorbar2onedose./2, '-go', 'LineWidth',2)
xlabel ('Time Post Treatment (Weeks)', 'FontSize', 16)
ylabel ('Fraction of Cells Resistant', 'FontSize', 16)
set(gca,'LineWidth',1.5,'FontSize',14)
errorbar(time, fres2untall, errorbar2untreated./2,'ko', 'LineWidth',2)
title('Two Population Model Effect of Pulse Treatments ','FontSize', 16)
legend('Twice Treated Resistant Fraction','Single Treatment Resistant Fraction', 'Untreated Resistant Fraction')

figure(5)
hold off
errorbar(time(1:10), fres2twodoseall(7:16), errorbar2twodose(7:16)./2, '-ro', 'LineWidth', 2)
hold on
errorbar(time(1:10),fres2onedoseall(1:10), errorbar2onedose(1:10)./2, '-go', 'LineWidth',2)
xlabel ('Time Post Treatment (Weeks)', 'FontSize', 16)
ylabel ('Fraction of Cells Resistant', 'FontSize', 16)
set(gca,'LineWidth',1.5,'FontSize',14)
errorbar(time(1:10), fres2untall(1:10), errorbar2untreated(1:10)./2,'ko', 'LineWidth',2)
title('Two Population Model Response to First and Second Treatment ','FontSize', 16)
legend('Resposne to Second Treatment','Response to First Treatment', 'Untreated Resistant Fraction')

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

%% Find Naive Overall Fractions
 % loads 3 column vector of dose, viability, and cohort number
 
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
v_model2allweeksnormednaive = model2popallweeksnormednaive( dose, coeffs2, beta2naive, Vmaxnaiveavg)

%% Perform Boot Strapping Estimates of Error in Parameters two population model
% Need to make this function and its naive counterpart
[lowerlim2boot, upperlim2boot] = BSerrorinparams2(residuals2, beta2newall, vddata);
[lowerlim2bootnaive, upperlim2bootnaive] = BSerrorinparams2naive(residuals2naive, dosenc, beta2naive, naive_vddata, beta2newall, Vmaxnaiveavg);

lowerlim2boot = lowerlim2boot
upperlim2boot = upperlim2boot
errorbarlength2boot = upperlim2boot-lowerlim2boot
lowerlim2bootnaive = lowerlim2bootnaive
upperlim2bootnaive = upperlim2bootnaive
errorbarlength2bootnaive = upperlim2bootnaive - lowerlim2bootnaive

%% Plot of resistant and sensitive fraction estimates over time

beta2unt = load('beta2unt.mat');
beta2unt = struct2cell(beta2unt);
beta2unt = cell2mat(beta2unt);

errorbarlength2bootuntreated = load('errorbarlength2bootuntreated.mat');
errorbarlength2bootuntreated = struct2cell(errorbarlength2bootuntreated);
errorbarlength2bootuntreated = cell2mat(errorbarlength2bootuntreated);
fres2 = zeros([15,1]);
fsens2 = beta2newall(5:19);

for i = 1:15
    fres2(i) = 1-fsens2(i);
end
time = zeros([1 15]);
time(1, :) = 1:15;

fres2unt = zeros([15,1]);
fsens2unt = beta2unt;

for i = 1:15
    fres2unt(i) = 1-fsens2unt(i);
end
time = zeros([1 15]);
time(1, :) = 1:15;

% Plot of the resistant and sensitive fractions
figure(6)
hold off
errorbar(time, fres2, errorbarlength2boot(5:19)./2, 'bo-', 'LineWidth', 1.5)
hold on
%errorbar(time, fsens2, errorbarlength2boot(5:19)./2, 'go-', 'LineWidth', 1.5)
errorbar(0, fres0, errorbarlength2bootnaive./2, 'bo', 'LineWidth', 1.5)
errorbar(time, fres2unt, errorbarlength2bootuntreated(5:19)./2, 'ro-', 'LineWidth', 2)
errorbar(0, fsens0, errorbarlength2bootnaive./2, 'ro', 'LineWidth', 1.5)
ylim([ 0 1])
xlabel('Time (Weeks Post Treatment)', 'FontSize', 14)
ylabel('Fraction of Cells', 'FontSize',14)
title ('Fraction of Resistant Cells Two Population Model', 'FontSize', 14)
set(gca,'LineWidth',1.5,'FontSize',14)
legend ('Fraction of Treated Cells Resistant', 'Fraction of Naive Cells Resistant', 'Fraction of Untreated Cells Resistant')
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

