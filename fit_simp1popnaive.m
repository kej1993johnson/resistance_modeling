function [unweighted_diffnew] = fit_simp1popnaive( beta1new, dose, var, Vmaxnaiveavg)
%This function is called by lsqnonlin.
% beta is a vector that contains the coefficeints of the equation. dose and
% var are the input data sets that were passed to lswqnonlin
slope1 = beta1new(1); % slope
cen1 = beta1new(2); % MD50 bulk
n= length(dose);


 
unweighted_diffnew = (Vmaxnaiveavg./( 1 + exp(slope1.*(dose - cen1)))) - var;



 
end