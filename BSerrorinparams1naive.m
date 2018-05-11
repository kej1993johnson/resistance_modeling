function [lowerlim1bootnaive, upperlim1bootnaive] = BSerrorinparams1naive(residualsLD50snaive_8, dosenc, betaLD50naive_8, cohort_8_naive, Vmaxnaiveavgcohort);

v_model1_naive=model1pop( dosenc, betaLD50naive_8, Vmaxnaiveavgcohort );
 options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
nboot = 500;
[~, bootIndices] = bootstrp(nboot, [], residualsLD50snaive_8); % randomly generates indices
bootResiduals = residualsLD50snaive_8(bootIndices); % uses indices to sample from residuals with replacement
varBoot = repmat(v_model1_naive,1,nboot) + bootResiduals; % creates simulated data set
% build up the bootstrap data sets to the single population function
betaBoot = zeros (nboot,2);
cohort_8_boot = cohort_8_naive;

for i = 1:nboot
    cohort_8_boot(:,3) = varBoot(:,i);
    paramslbn = zeros([2 1]);
dose0ind = cohort_8_boot(:,2) == 0;
 Vmaxcohortnaivedata = cohort_8_boot(dose0ind,:);
 Vmaxnaiveavgcohort_boot = mean(Vmaxcohortnaivedata(:,3));
paramsubn = Inf([2 1]);

params0n = [ .1; 30];
% outputs betaLD50 with first 12 paramters slope values, second 12 are
% center values in week order
betaBoot(i,:) = lsqnonlin(@fit_simp1popnaive,...
    params0n,...
    paramslbn,...
    paramsubn,...
    options,...
    dosenc,...
    varBoot(:,i),...
    Vmaxnaiveavgcohort_boot);
    
   
end


 bootCI = prctile(betaBoot, [2.5 97.5]);
 lowerlim1bootnaive = bootCI(1,:);
 upperlim1bootnaive = bootCI(2,:);
end