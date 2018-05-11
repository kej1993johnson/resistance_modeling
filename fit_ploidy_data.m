 clear all, close all, clc

vddata = load ('ploidy_3wpt_data.m'); % load concatenated data set (long Mx2 matrix
% all weeks added vertically, first column is week number
dose = vddata (:,2);
var = vddata (:,3);
WPT = vddata(:,4);
n = length(dose); % finds number of data points total
ploidy_num = vddata(:,1); % note this is allocated as 1 2 3 4 5 which corresponds to untreated, wp, 2n, 4n, 8n
[ Vmaxbyweek, Vmaxweekavg, ninweek, wknum, Vmaxall] = findVmaxandsize(vddata);

%% Two Population Model
options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
%define bounds for all parameters
            % m1 cen1 m2 cen 2 wk1sens wk1res etc...
paramslb = zeros([9 1]);
paramsub = [ Inf; Inf; Inf; Inf; 1; 1; 1; 1; 1];

% sets up initial values of beta3new
params0 = [ .1; 17; .04; 35; .3; .3; .3; .3; .3];

%  res2sq=ones([n 1]);
[beta2ploidy, resnormp, residualsp] = lsqnonlin(@fit_ploidy8N_2pop,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    dose,...
    var,...
    wknum,...
    Vmaxall);
coeffs2ploidy = beta2ploidy(1:4);
%%
for i = 1:5
fsens(i) = beta2ploidy(4+(i));
end

for i = 1:5
fres(i)= 1-fsens(i);
end

time = 1:1:5;

figure(11)
hold off
hold on
plot(time(1), fsens(1), 'yo', 'LineWidth', 2.5)
plot(time(1), fres(1), 'bo', 'LineWidth', 2.5)
plot(time(2), fsens(2), 'yo', 'LineWidth', 2.5)
plot(time(2), fres(2), 'bo', 'LineWidth', 2.5)
plot(time(3:5), fsens(3:5), '-yo', 'LineWidth', 2.5)
plot(time(3:5), fres(3:5), '-bo', 'LineWidth', 2.5)
ylim([0 1])
xlabel ('Ploidy (Untreated, Whole Population 3WPT, 2N, 4N, 8N', 'FontSize', 16)
ylabel ('Fraction of Cells', 'FontSize', 16)
title (' Two Population Model on Ploidy Data', 'FontSize', 16)
legend ('Sensitive Fraction LD50 = 38.8 uM','Resistant Fraction LD50 = 81.7 uM')
set(gca,'LineWidth',1.5,'FontSize',16)


 %% Three Population Model Fit 
 options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
%define bounds for all parameters
            % m1 cen1 m2 cen 2 wk1sens wk1res etc...
paramslb = zeros([16 1]);
paramsub = [ Inf; Inf; Inf; Inf; Inf; Inf; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1];

% sets up initial values of beta3new
params0 = [ .1; 17; .04; 35; .04; 65; .3; .3; .3; .3; .3; .3; .3; .3; .3; .3];

%  res2sq=ones([n 1]);
[beta3ploidy, resnormp, residualsp] = lsqnonlin(@fit_ploidy8N_3pop,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    dose,...
    var,...
    wknum,...
    Vmaxall);
coeffs3ploidy = beta3ploidy(1:6);


for i = 1:5
fsens(i) = beta3ploidy(6+(2.*i-1));
end

for i = 1:5
    fres(i)= beta3ploidy(6+(2.*i));
end

for i = 1:5
    ftol(i) = 1- fres(i)-fsens(i);
end
time = 1:1:5;
%%
%%% Note need to plot using bar graphs! Or in r...
fres_new = ftol;
ftol_new = fres;

figure(11)
hold off
hold on
plot(time(1), fsens(1), 'yo', 'LineWidth', 2.5)
plot(time(1), ftol_new(1), 'go', 'LineWidth', 2.5)
plot(time(1), fres_new(1), 'bo', 'LineWidth', 2.5)
plot(time(2), fsens(2), 'yo', 'LineWidth', 2.5)
plot(time(2), ftol_new(2), 'go', 'LineWidth', 2.5)
plot(time(2), fres_new(2), 'bo', 'LineWidth', 2.5)
plot(time(3:5), fsens(3:5), '-yo', 'LineWidth', 2.5)
plot(time(3:5), ftol_new(3:5), '-go', 'LineWidth', 2.5)
plot(time(3:5), fres_new(3:5), '-bo', 'LineWidth', 2.5)
ylim([0 1])
xlabel ('Ploidy (Untreated, Whole Population, 2N, 4N, 8N', 'FontSize', 16)
ylabel ('Fraction of Cells', 'FontSize', 16)
title (' Three Population Model on Ploidy Data', 'FontSize', 16)
legend ('Sensitive Fraction LD50 = 9.9 uM', 'Tolerant Fraction LD50 = 51.6 uM ',' Resistant Fraction LD50 = 89.7 uM')
set(gca,'LineWidth',1.5,'FontSize',16)

