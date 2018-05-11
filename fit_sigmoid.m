function [unweighted_diffnew] = fit_sigmoid( beta1new, dose, var, Vmaxall )
%This function is called by lsqnonlin.
% beta is a vector that contains the coefficeints of the equation. dose and
% var are the input data sets that were passed to lswqnonlin
slope1 = beta1new(1); % slope
cen1 = beta1new(2); % MD50 bulk


 
unweighted_diffnew = (Vmaxall./( 1 + exp(slope1.*(dose - cen1)))) - var;

%weighted_diffnew = unweighted_diffnew.* new_vddata(:,5);



 
end