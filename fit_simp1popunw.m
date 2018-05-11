function [unweighted_diffnew] = fit_simp1popunw( beta1new, dose, var, wknum, Vmaxall )
%This function is called by lsqnonlin.
% beta is a vector that contains the coefficeints of the equation. dose and
% var are the input data sets that were passed to lswqnonlin
slope1 = beta1new(1); % slope
cen1 = beta1new(2); % MD50 bulk
n= length(dose);

nsize = wknum(:,2);


n= length(dose);

 
unweighted_diffnew = (Vmaxall./( 1 + exp(slope1.*(dose - cen1)))) - var;

%weighted_diffnew = unweighted_diffnew.* new_vddata(:,5);



 
end