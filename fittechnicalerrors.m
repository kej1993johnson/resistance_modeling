function [kstd, Wstd, std_dev_model_fun] = fittechnicalerrors(dose)

% Fits the same exponential function except instead of using residuals uses
% error in technical replicates at each dose. Since doses are the same,
% binning as in above isn't required and function is much simpler. One
% issue is that this may not be extrapolated to be the proper error as a
% function of dose
n = length(dose);

dose_res_avg_all = load('standard_dev_in_tech_error.m');

doseavg = dose_res_avg_all(:,1);
res1avg = dose_res_avg_all(:,2);

k1guess = [0.01; .03; 15];

[k1] = lsqnonlin(@fitabsres1, k1guess,[0; 0; 0],[ Inf; Inf; Inf],[],doseavg,res1avg);

kstd = k1;
%exp(slope1.*(dose - cen1)))
X = (.1:.1:250);
std_dev_model_fun = (k1(1)*(X+k1(3)).*exp(-k1(2)*X));
Wstd = 1./(std_dev_model_fun);


end
