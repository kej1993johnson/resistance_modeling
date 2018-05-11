clear all, close all, clc

vddata = load ('allweeksdata11_26_14.m'); % load concatenated data set (long Mx2 matrix
% all weeks added vertically, first column is week number

%variance = var(vddata)
week = vddata(:,1);
dose = vddata (:,2);
var = vddata (:,3);
cohort_number = vddata(:,4);

n = length(dose); % finds number of data points total

[ Vmaxbyweek, Vmaxweekavg, ninweek, wknum, Vmaxall] = findVmaxandsize(vddata);
nsize = wknum(:,2);

%% Finds weighted residual functions

[k1,k2] = fitabsoluteresiduals14w(dose, var, wknum, Vmaxall);

[kstd, Wstd, std_dev_model_fun] = fittechnicalerrors(dose);

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
paramslb = zeros([18 1]);
paramsub = [ Inf; Inf; Inf; Inf; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1];

% sets up initial values of beta2new
params0 = [ .1; 17; .04; 50; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5];
%  res2sq=ones([n 1]);
[beta2new, resnorm2, residuals2] = lsqnonlin(@fit_simp2popabsresdensnormed14w,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    dose,...
    var,...
    k2,...
    wknum,...
    Vmaxall);
%%
v_model2allweeksnormed = model2popallweeksnormed( dose, beta2new, Vmaxbyweek, nsize);
actualresiduals2 = var- v_model2allweeksnormed;
sigma2 = std(actualresiduals2)
error2 = sum(abs(actualresiduals2))./n
chi_squared2 = sum((actualresiduals2.^2)./v_model2allweeksnormed)
RSS2 = sum(actualresiduals2.^2);
p2 = 16;
AIC2 = n*log(RSS2./n) + 2*p2
BIC2 = -n*log(RSS2./n) + p2*log(n)

coeffs2 = beta2new(1:4)

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
%[lowerlim2, upperlim2, errorbarlength2, beta2newsim ] = MCerrorinparams2(dose, beta2new, Vmaxbyweek, nsize, sigma2, k2, kstd, wknum);

%%
lowerboundsparams2 = lowerlim2(1:4)
upperboundsparams2 = upperlim2(1:4)
%%
fres2 = zeros( [14,1]);
fsens2 = beta2new(5:18);

for i = 1:14
    fres2(i) = 1-fsens2(i);
end
time = zeros([1 14]);
time(1, :) = 1:14;

% Find error bars

% Plot of the resistant and sensitive fractions
figure(2)
hold off
plot(time, fres2, 'b', 'LineWidth', 4)
%errorbar(time, fres2, errorbarlength2(5:16)./2, 'b', 'LineWidth', 3)
hold on
plot(time, fsens2, 'g', 'LineWidth',4)
%errorbar(time, fsens2, errorbarlength2(5:16)./2, 'g', 'LineWidth',3)
legend('resistant (LD50 = 23.3 +/- 1.95 uM)', 'sensitive (LD50 = 81.1 +/- 3.65 uM)')
xlim([ 1 14])
ylim([ 0 1])
xlabel('Weeks Post Treatment', 'FontSize', 24)
ylabel('Fraction of Cells', 'FontSize',24)
%title ('Fraction of Resistant Cells and Sensitive Cells vs. Time Post Treatment', 'FontSize', 14)
hold off

%%

 % calculates AIC value for new model. Again, not sure about the p-value
 % here since we're technically fitting this data to 28 betas
 % (slope1, xcen1, slope2, xcen2, +12 weeks with two fractions in each week)


% Function creates vector of predicted viabilities based on model for each
% week
% First is for the unnormalized model
v_model2_3 = model2pop16p( X, beta2new, Vmaxbyweek, 3);
v_model2_6 = model2pop16p( X, beta2new, Vmaxbyweek, 6);
v_model2_10 = model2pop16p( X, beta2new, Vmaxbyweek, 10);

%%
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

v_model2allweeksnormedfinal = model2popallweeksnormed( dose, params2new, Vmaxbyweek, nsize);
figure(4)
hold off
plot(dose,v_model2allweeksnormedfinal, 'ko', 'LineWidth', 4)
hold on
plot(dose, var, '.')
xlim ([0 250])
xlabel('dose (uM)','FontSize',18)
ylabel('viability','FontSize',18)
title ('Two Population AIC Weighted Final Model','FontSize',18)


