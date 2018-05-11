function [lowerlim1boot, upperlim1boot] = BSerrorinparams1naivebleh(residuals1, dose, beta1_8, Vmaxweekavg, cohort_8)
% Next try bootstrapping method...
v_model1=model1pop( dose, beta1_8, Vmaxweekavg );
 options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
nboot = 500;
[~, bootIndices] = bootstrp(nboot, [], residuals1); % randomly generates indices
bootResiduals = residuals1(bootIndices); % uses indices to sample from residuals with replacement
varBoot = repmat(v_model1,1,nboot) + bootResiduals; % creates simulated data set
% build up the bootstrap data sets to the single population function
betaBoot = zeros (nboot,2);
cohort_8_boot = cohort_8;
for i = 1:nboot
    paramslbn = zeros([2 1]);
paramsubn = Inf([2 1]);
dose0ind = cohort_8_boot(:,2) == 0;
 Vmaxcohortnaivedata = cohort_8_boot(dose0ind,:);
 Vmaxnaiveavgcohort_boot = mean(Vmaxcohortnaivedata(:,3));
params0n = [ .1; 30];
% outputs betaLD50 with first 12 paramters slope values, second 12 are
% center values in week order
[ Vmaxbyweek, Vmaxweekavg_boot, ninweek, wknum_boot, Vmaxall_boot] = findVmaxandsizetest15(cohort_8_boot);
betaBoot(i,:) = lsqnonlin(@fit_simp1popnaive,...
    params0n,...
    paramslbn,...
    paramsubn,...
    options,...
    dose,...
    varBoot(:,1),...
    Vmaxnaiveavgcohort_boot);
    
   
end


 bootCI = prctile(betaBoot, [2.5 97.5]);
 lowerlim1boot = bootCI(1,:);
 upperlim1boot = bootCI(2,:);
end
