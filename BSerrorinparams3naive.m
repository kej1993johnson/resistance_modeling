function [lowerlim3bootnaive, upperlim3bootnaive] = BSerrorinparams3naive(residuals3naive, dosenc, beta3naive, naive_vddata, beta3new, Vmaxnaiveavg)
coeffs3 = beta3new(1:6);
v_model3_naive = model3popallweeksnormednaive( dosenc, coeffs3, beta3naive);
 options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
nboot = 500;
[~, bootIndices] = bootstrp(nboot, [], residuals3naive); % randomly generates indices
bootResiduals = residuals3naive(bootIndices); % uses indices to sample from residuals with replacement
varBoot = repmat(v_model3_naive,1,nboot) + bootResiduals; % creates simulated data set
% build up the bootstrap data sets to the single population function

vddata_boot = naive_vddata;

for i = 1:nboot
    vddata_boot(:,3) = varBoot(:,i);
    paramslbn = zeros([2 1]);
dose0ind = vddata_boot(:,2) == 0;
 Vmaxnaivedata = vddata_boot(dose0ind,:);
 Vmaxnaiveavg_boot = mean(Vmaxnaivedata(:,3));
 
 params0 = .5;
paramslb = 0;
paramsub = 1;

betabootnaive(i) = lsqnonlin(@fit_simp3pop6weeksnaive,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    dosenc,...
    varBoot(:,i),...
    coeffs3);
 
% paramsubn = Inf([2 1]);
% 
% params0n = [ .1; 30];
% % outputs betaLD50 with first 12 paramters slope values, second 12 are
% % center values in week order
% betaBoot(i,:) = lsqnonlin(@fit_simp1popnaive,...
%     params0n,...
%     paramslbn,...
%     paramsubn,...
%     options,...
%     dosenc,...
%     varBoot(:,i),...
%     Vmaxnaiveavgcohort_boot);
    
   
end


 lowerlim2bootnaive = prctile(betabootnaive', 2.5);
 upperlim2bootnaive = prctile(betabootnaive', 97.5);
end
