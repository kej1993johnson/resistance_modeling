function [ lowerlim, upperlim, errorbarlength, beta2newsim ] = MCerrorinparams2( vddata, v_model2allweeksnormed, beta2, Vmaxbyweek, kstd, wknum, Vmaxall)

% MCerrorinparamsLD50( vddata, bulk_model_LD50, betaLD50, kstd, wknum, Vmaxall)
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

% v_model2allweeks= Vmaxbyweek(1).*(beta2new(5)./( 1 + exp(beta2new(1).*(dose(1:sum(nsize(1))) - beta2new(2))))) + (beta2new(6)./(1 + exp(beta2new(3).*(dose(1:sum(nsize(1))) - beta2new(4)))));
% for i = 2: 12
% v_model2 = Vmaxbyweek(i).*(beta2new(4+i)./( 1 + exp(beta2new(1).*(dose(sum(nsize(1:(i-1)))+1:sum(nsize(1:i))) - beta2new(2))))) + (beta2new(4+i)./(1 + exp(beta2new(3).*(dose(sum(nsize(1:(i-1)))+1:sum(nsize(1:i))) - beta2new(4)))));
% v_model2allweeks = vertcat(v_model2allweeks, v_model2);
% end
nsize = wknum(:,2);
dose = vddata(:,2);
% 
% v_model2allweeksnormed= Vmaxbyweek(1).*(beta2new(5)./( 1 + exp(beta2new(1).*(dose(1:sum(nsize(1))) - beta2new(2))))) + ((1-beta2new(5))./(1 + exp(beta2new(3).*(dose(1:sum(nsize(1))) - beta2new(4)))));
% 
%  for i = 2: 12
% %v_model2 = Vmaxbyweek(i).*(beta2new(3+2*i)./( 1 + exp(beta2new(1).*(dose(sum(nsize(1:(i-1)))+1:sum(nsize(1:i))) - beta2new(2))))) + (beta2new(4+2*i)./(1 + exp(beta2new(3).*(dose(sum(nsize(1:(i-1)))+1:sum(nsize(1:i))) - beta2new(4)))));     
% v_model2 = Vmaxbyweek(i).*(beta2new(4+i)./( 1 + exp(beta2new(1).*(dose(sum(nsize(1:(i-1)))+1:sum(nsize(1:i))) - beta2new(2))))) + (1-beta2new(4+i))./(1 + exp(beta2new(3).*(dose(sum(nsize(1:(i-1)))+1:sum(nsize(1:i))) - beta2new(4))));
% v_model2allweeksnormed = vertcat(v_model2allweeksnormed, v_model2);
%  end

% Now set this equal to the new data set
datasetnew = v_model2allweeksnormed;

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
% Now find paramters for the simulated data set

options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);

paramslb = zeros([19 1]);
paramsub = [ Inf; Inf; Inf; Inf; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1;];

% sets up initial values of beta2new
params0 = [ .1; 17; .04; 50; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5];
beta2newsim(:,i) = lsqnonlin(@fit_simp2popunw15,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    dose,...
    datasetsim,...
    wknum,...
    Vmaxall);


end


beta2newsimrot = beta2newsim';
stdall = std(beta2newsimrot);
%stdallcorr = stdall'./(sqrt(nruns)); % this is s(x)/sqrt(n) = sigma(sample mean)
% Consider changing this because ndose is very large, might actually depend
% on nruns where this is what follows the Gaussian distribution around the
% "sample mean" which is the parameter


stdallcorr = stdall';

% now turn this into a confidence interval over beta2new where beta2new(i)
% is the mean and the 95 percent confidence interval is the upper and lower
% limit of the mean +/- 2* std(beta2newsim(i)

lambda = 1.96;

lowerlim = zeros([19,1]);
upperlim = zeros([19,1]);
errorbarlength = zeros([19,1]);
for i = 1:19
lowerlim(i, 1) = beta2(i) - lambda.*stdallcorr(i);
upperlim(i,1) = beta2(i) + lambda.*stdallcorr(i);
errorbarlength(i,1) = upperlim(i,1) - lowerlim(i,1);
end
end

