function [lowerlim2, upperlim2] = BSerrorinparams2bleh(residuals2, beta2, cohort_8)
% Next try bootstrapping method...
dose = cohort_8(:,2);

[ Vmaxbyweek, Vmaxweekavg, ninweek, wknum, Vmaxall] = findVmaxandsizetest15(cohort_8);
nsize = wknum(:,2);
v_model2=model2popallweeksnormed( dose, beta2, Vmaxbyweek, nsize);
nboot = 10;
[~, bootIndices] = bootstrp(nboot, [], residuals2); % randomly generates indices
bootResiduals = residuals2(bootIndices); % uses indices to sample from residuals with replacement
varBoot = repmat(v_model2,1,nboot) + bootResiduals; % creates simulated data set
% build up the bootstrap data sets to the single population function
betaboot = zeros (nboot,2);
cohort_8_boot = cohort_8;
options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
%define bounds for all parameters
            % m1 cen1 m2 cen 2 wk1sens wk1res etc...
paramslb = zeros([19 1]);
paramsub = [ Inf; Inf; Inf; Inf; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1;];

% sets up initial values of beta2new
params0 = [ .1; 17; .04; 50; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5];
for i = 1:nboot
    cohort_8_boot(:,3) = varBoot(:,i);
    [ Vmaxbyweek, Vmaxweekavg_boot, ninweek_boot, wknum_boot, Vmaxall_boot] = findVmaxandsizetest15(cohort_8_boot);
    betaboot(:,i)= lsqnonlin(@fit_simp2popunw15,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    dose,...
    varBoot(:,i),...
    wknum_boot,...
    Vmaxall_boot);
     
end


 bootCI = prctile(betaboot, [2.5 97.5]);
 lowerlim2 = bootCI(1,:);
 upperlim2 = bootCI(2,:);
end
