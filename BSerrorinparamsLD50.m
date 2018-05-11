function [lowerlim1boot, upperlim1boot] = BSerrorinparamsLD50(residuals, vddata, betaLD50)
% Next try bootstrapping method...

dose = vddata(:,2);
wk = vddata(:,1);
[ Vmaxbyweek, Vmaxweekavg, ninweek, wknum, Vmaxall] = findVmaxandsize(vddata);
nsize = wknum(:,2);
v_model1 = model1popallweeksnormed( dose, betaLD50, Vmaxbyweek, nsize );

nboot = 500;
[~, bootIndices] = bootstrp(nboot, [], residuals); % randomly generates indices
bootResiduals = residuals(bootIndices); % uses indices to sample from residuals with replacement
varBoot = repmat(v_model1,1,nboot) + bootResiduals; % creates simulated data set
% build up the bootstrap data sets to the single population function
% betaBoot = zeros(nboot,30);
vddata_boot = vddata;
for i = 1:nboot
    vddata_boot(:,3) = varBoot(:,i);
    
    [ Vmaxbyweek_boot, Vmaxweekavg_boot, ninweek, wknum_boot, Vmaxall_boot] = findVmaxandsize(vddata_boot);
if length(Vmaxbyweek) == 15 
paramslb = zeros([30 1]);
paramsub = Inf([30 1]);

options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);

params0 = [ .1; .1; .1; .1; .1;.1;.1; .1; .1; .1; .1; .1; .1; .1; .1; 30; 30; 30; 30; 30; 30; 30; 30; 30; 30; 30; 30; 30; 30; 30];
% outputs betaLD50 with first 12 paramters slope values, second 12 are
% center values in week order

betaboot(i,:) = lsqnonlin(@fit_simpLD5015,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    dose,...
    varBoot(:,i),...
    wknum_boot,...
    Vmaxall_boot);
end 

if length(Vmaxbyweek) == 8 
paramslb = zeros([16 1]);
paramsub = Inf([16 1]);

options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);

params0 = [ .1; .1; .1; .1; .1;.1;.1; .1; 30; 30; 30; 30; 30; 30; 30; 30];
% outputs betaLD50 with first 12 paramters slope values, second 12 are
% center values in week order

betaboot(i,:) = lsqnonlin(@fit_simpLD508,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    dose,...
    varBoot(:,i),...
    wknum_boot,...
    Vmaxall_boot);
end 
end


 bootCI = prctile(betaboot, [2.5 97.5]);
 lowerlim1boot = bootCI(1,:);
 upperlim1boot = bootCI(2,:);
 end