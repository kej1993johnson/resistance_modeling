function[lowerlim3, upperlim3] = BSerrorinparams3(residuals3, beta3new, vddata)
dose = vddata(:,2);


[ Vmaxbyweek, Vmaxweekavg, ninweek, wknum, Vmaxall] = findVmaxandsize(vddata);
nsize = wknum(:,2);
v_model3 = model3popallweeksnormed( dose, beta3new, Vmaxbyweek, nsize);
nboot = 500;
[~, bootIndices] = bootstrp(nboot, [], residuals3); % randomly generates indices
bootResiduals = residuals3(bootIndices); % uses indices to sample from residuals with replacement
varBoot = repmat(v_model3,1,nboot) + bootResiduals; % creates simulated data set
% build up the bootstrap data sets to the single population function
vddata_boot = vddata;
options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
%define bounds for all parameters
            % m1 cen1 m2 cen 2 wk1sens wk1res etc...
paramslb = zeros([22 1]);
paramsub = [ Inf; Inf; Inf; Inf; Inf; Inf; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1];
%%
% sets up initial values of beta2new
params0 = [ .1; 17; .04; 50; .3; 80; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5];
for i = 1:nboot
    vddata_boot(:,3) = varBoot(:,i);
    [ Vmaxbyweek_boot, Vmaxweekavg_boot, ninweek_boot, wknum_boot, Vmaxall_boot] = findVmaxandsize(vddata_boot);
    betaboot(:,i)= lsqnonlin(@fit_simp3pop8weeks,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    dose,...
    varBoot(:,i),...
    wknum_boot,...
    Vmaxall_boot);
     
end


 bootCI = prctile(betaboot', [2.5 97.5]);
 lowerlim3 = bootCI(1,:);
 upperlim3 = bootCI(2,:);
end