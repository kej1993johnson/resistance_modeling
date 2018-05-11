% This script is for the second treatment at 2 weeks and does analysis out
% to 8 weeks 

clear all, close all, clc

%vddata = load ('secondtreatment.m'); % only dosed once ones
vddata = load ('second_treatment_at_2_wks_data.m'); % dosed twice
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

%% Two Dose Two population model

coeffs2 = load('coeffs2.mat');
coeffs2 = struct2cell(coeffs2);
coeffs2 = cell2mat(coeffs2);
options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
%define bounds for all parameters
            % m1 cen1 m2 cen 2 wk1sens wk1res etc...
paramslb = zeros([6 1]);
paramsub = ones([6 1]);

% sets up initial values of beta2new
params0 = 0.3.* ones([6 1]);
%  res2sq=ones([n 1]);
[beta2twodose, resnorm2, residuals2] = lsqnonlin(@fit_simp2popunw6_untreated,...
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
beta2onedose = coeffs2;
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


fres2twodose = zeros([6,1]);
fsens2twodose = beta2twodose;

for i = 1:6
    fres2twodose(i) = 1-fsens2twodose(i);
end

fres2twodoseall = vertcat(fres0, fres2twodose);
fsens2twodoseall = vertcat(fsens0, fsens2twodose);

errorbarlength2bootuntreated = load('errorbarlength2bootuntreated.mat');
errorbarlength2bootuntreated = struct2cell(errorbarlength2bootuntreated);
errorbarlength2bootuntreated = cell2mat(errorbarlength2bootuntreated);
errorbar2untreated = vertcat( errorbarlength2bootnaive, errorbarlength2bootuntreated(5:19)');
%errorbar2untreated is 16 weeks of fraction estimates

%errorbar2onedose = vertcat(errorbarlength2bootnaive, errorbarlength2onedose(5:10)');
errorbar2twodose = vertcat(errorbarlength2bootnaive, errorbarlength2twodose(5:10)');

%%
time = 0:1:6;
figure (4)
hold off
% Plot of LD50s over time
%plot(time, betaLD50(16:30), 'o', 'LineWidth', 4)
errorbar(time, fres2twodoseall, errorbar2twodose./2, '-bo', 'LineWidth', 2)
hold on
errorbar(time, fsens2twodoseall, errorbar2twodose./2, '-yo', 'LineWidth', 2)
xlabel ('Time Post Treatment (Weeks)', 'FontSize', 16)
ylabel ('Fraction of Cells Resistant', 'FontSize', 16)
set(gca,'LineWidth',1.5,'FontSize',14)
%errorbar(time, fres2untall, errorbar2untreated./2,'ko', 'LineWidth',2)
title('Two Population Model Effect of Pulse Treatments ','FontSize', 16)
legend('Twice Treated Resistant Fraction','Single Treatment Resistant Fraction', 'Untreated Resistant Fraction')

