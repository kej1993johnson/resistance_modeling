function [lowerlim1boot, upperlim1boot, AIC, MSE] = BSerrorinparams1(residuals1, dose, beta1_8, Vmaxweekavg, cohort_8)
% Next try bootstrapping method...
v_model1=model1pop( dose, beta1_8, Vmaxweekavg );
nboot = 500;
[~, bootIndices] = bootstrp(nboot, [], residuals1); % randomly generates indices
bootResiduals = residuals1(bootIndices); % uses indices to sample from residuals with replacement
varBoot = repmat(v_model1,1,nboot) + bootResiduals; % creates simulated data set
% build up the bootstrap data sets to the single population function
betaBoot = zeros (nboot,2);
cohort_8_boot = cohort_8;
initials1new = [ .05; 30];
n = length(dose);
p1 = 2;
for i = 1:nboot
    cohort_8_boot(:,3) = varBoot(:,i);
    [ Vmaxbyweek, Vmaxweekavg_boot, ninweek, wknum_boot, Vmaxall_boot] = findVmaxandsizetest15(cohort_8_boot);
    [betaBoot(i,:), residuals] = lsqnonlin(@fit_simp1popunw, initials1new,[0; 0],[ Inf; Inf],[],dose,varBoot(:,i), wknum_boot, Vmaxall_boot);
    RSS = sum(residuals.^2);
    AIC(i) = n*log(RSS./n) + 2*p1;
    v_model=model1pop( dose, betaBoot(i,:), Vmaxweekavg_boot);
    chi_squared = sum((residuals.^2)./v_model);
    DOF = n-p1;
    MSE(i) = chi_squared./ DOF;
end


 bootCI = prctile(betaBoot, [2.5 97.5]);
 lowerlim1boot = bootCI(1,:);
 upperlim1boot = bootCI(2,:);
end
