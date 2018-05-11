function difforignew = fit_ploidy8N_2pop( beta2new, dose, var, wknum, Vmaxall )
%This function is called by lsqnonlin.
% beta is a vector that contains the coefficeints of the equation. dose and
% var are the input data sets that were passed to lswqnonlin
slope1 = beta2new(1); % slope
cen1 = beta2new(2); % MD50 bulk
slope2 = beta2new(3);
cen2 = beta2new(4);
f11 = beta2new(5);
f21 = beta2new(6);
f31 = beta2new(7);
f41 = beta2new(8);
f51 = beta2new(9);

n= length(dose);

nsize = wknum(:,2);

fvec1 = zeros([n 1]);
fvec1(1:sum(nsize(1))) = f11;
fvec1((sum(nsize(1))+1:sum(nsize(1:2)))) = f21;
fvec1((sum(nsize(1:2))+1:sum(nsize(1:3)))) = f31;
fvec1((sum(nsize(1:3))+1:sum(nsize(1:4)))) = f41;
fvec1((sum(nsize(1:4))+1:sum(nsize(1:5)))) = f51;




 difforignew = ((((fvec1)./( 1 + exp(slope1.*(dose - cen1)))) + ((1-fvec1)./(1 + exp(slope2.*(dose - cen2))))).*Vmaxall) - var;