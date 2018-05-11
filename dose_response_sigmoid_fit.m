 %% Find LD50 from Single Dose Response Curve
 close all, clear all, clc
 
% dose_resp_data_ADR = load ('dose_resp_ADR_7.m') 
% dose_resp_data_WT = load('dose_resp_MCF7_7.m')
%  
dose_resp_data_ADR = load ('dose_resp_ADR_24h.m'); 
dose_resp_data_WT = load('dose_resp_MCF7_24h.m');
dose0ind = dose_resp_data_ADR(:,1) == 0;
 VmaxallADR = dose_resp_data_ADR(dose0ind,:);
 VmaxallavgADR = mean(VmaxallADR(:,2));
 
 dose = dose_resp_data_ADR(:,1);
 viability = dose_resp_data_ADR(:,2);
 initials1new = [ .05; 30];
  options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
 [betaADR, resnorm, residuals] = lsqnonlin(@fit_sigmoid, initials1new,[0; 0],[ Inf; Inf],[],dose,viability, VmaxallavgADR);
X = (0:.1:400);
X = X';
 v_modelplotADR = model1pop(X, betaADR, VmaxallavgADR);

 %%
 dose0ind = dose_resp_data_WT(:,1) == 0;
 VmaxallWT = dose_resp_data_WT(dose0ind,:);
 VmaxallavgWT = mean(VmaxallWT(:,2));
 
 doseWT = dose_resp_data_WT(:,1);
 viabilityWT = dose_resp_data_WT(:,2);
 initials1new = [ .05; 30];
  options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
 [betaWT, resnormWT, residualsWT] = lsqnonlin(@fit_sigmoid, initials1new,[0; 0],[ Inf; Inf],[],doseWT,viabilityWT, VmaxallavgWT);
X = (0:.1:400);
X = X';
 v_modelplotWT = model1pop(X, betaWT, VmaxallavgWT);
 
 %%
 
 figure(1) % graphs fit to all weeks in red
hold off
plot (dose, viability, 'bo')
hold on
plot (doseWT, viabilityWT, 'ro')
plot (X, v_modelplotADR, '-b', 'LineWidth', 3);
plot( X, v_modelplotWT, '-r', 'LineWidth', 3);
xlim([0,400]);
xlabel('Dose (uM)','FontSize',18)
ylabel('Viability','FontSize',18)
title ('MCF7-ADR and WT 24 hour Dox Tolerance Test','FontSize',18)
legend('ADR','WT', 'ADR Best Fit, LD50 = 270.2 uM', 'WT Best Fit, LD50 = 38.7 uM')
 
%%
