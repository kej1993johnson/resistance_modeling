function [ lowerlim, upperlim, errorbarlength, beta1sim ] = MCerrorinparamsnaive( vddata, v_model1, beta1, kstd)
% This function takes the "new data set" given by the model fit for each
% dose value, and creates n = 1000 data set (i)s which follow a Gaussian
% distribution with spread given by the standard deviation of the residuals
% of the original data set.  For each new data set, it calculates the
% paramaters that the lsqnonlin returns for that data set.

% Once we have the vector of 1000 of each parameter, we can find the
% standard deviation of each paramter in order to find the determine the
% 95% confidence interval of each paramter.  This will allow us to apply
% horizontal error bars at the centers of each single time course graph, as
% well as vertical error bars to the resistant and sensitive fractions over
% time

% First find the "new data set" given by the model

dose = vddata(:,2);

dose0ind = vddata(:,2) == 0;
 Vmaxcohortnaivedata = vddata(dose0ind,:);
 Vmaxnaiveavgcohort = mean(Vmaxcohortnaivedata(:,3));
% Now set this equal to the new data set
datasetnew = v_model1;

ndose = length(dose);
nruns = 10;
% first create noise vector to be added to new data set (should be length
% dose, spread given by sigma (standard deviation of resiudals)

% should noise be scaled by nweek for the fractional parameters? So we have
% more noise in these paramters since they are being found using only the
% data from that specific week


for i = 1:nruns

std_dev_model_fun = (kstd(1).*(dose+kstd(3)).*exp(-kstd(2)*dose));

noise = (std_dev_model_fun.*randn(ndose,1));
datasetsim = datasetnew + noise;
vddata_MC = vddata;
% Now find paramters for the simulated data set
vddata_MC(:,3) = datasetsim;
dose0ind = vddata_MC(:,2) == 0;
 Vmaxcohortnaivedata = vddata_MC(dose0ind,:);
 Vmaxnaiveavgcohort = mean(Vmaxcohortnaivedata(:,3));
%[ Vmaxbyweek_MC, Vmaxweekavg_MC, ninweek, wknum_MC, Vmaxall_MC] = findVmaxandsizetest15(vddata_MC);
options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);
paramslb = zeros([30 1]);
paramsub = Inf([30 1]);
initials1new = [ .05; 30]; % sets up initial values of beta1new

% call to lsqnonlin, which calls fit_simp1popnew which contains a
% difference vector in which the Vmax is replaced by different variables
% corresponding to the week. Here this parameter should be 1 every time
% (but its not)

beta1sim(:,i) = lsqnonlin(@fit_simp1popnaive, initials1new,[0; 0],[ Inf; Inf],options,dose,datasetsim, Vmaxnaiveavgcohort);


end


beta1simrot = beta1sim';
stdall = std(beta1simrot);
%stdallcorr = stdall'./(sqrt(nruns)); % this is s(x)/sqrt(n) = sigma(sample mean)
% Consider changing this because ndose is very large, might actually depend
% on nruns where this is what follows the Gaussian distribution around the
% "sample mean" which is the parameter


stdallcorr = stdall';

% now turn this into a confidence interval over beta2new where beta2new(i)
% is the mean and the 95 percent confidence interval is the upper and lower
% limit of the mean +/- 2* std(beta2newsim(i)

lambda = 1.96;

lowerlim = zeros([2,1]);
upperlim = zeros([2,1]);
errorbarlength = zeros([2,1]);
for i = 1:2
lowerlim(i, 1) = beta1(i) - lambda.*stdallcorr(i);
upperlim(i,1) = beta1(i) + lambda.*stdallcorr(i);
errorbarlength(i,1) = upperlim(i,1) - lowerlim(i,1);
end
end

