function difforignew = fit_simp3pop6weeks_untreated( beta3new, dose, var, wknum, Vmaxall,coeffs3 )
%This function is called by lsqnonlin.
% beta is a vector that contains the coefficeints of the equation. dose and
% var are the input data sets that were passed to lswqnonlin
slope1 = coeffs3(1); % slope
cen1 = coeffs3(2); % MD50 bulk
slope2 = coeffs3(3);
cen2 = coeffs3(4);
slope3 = coeffs3(5);
cen3 = coeffs3(6);
f11 = beta3new(1);
f12 = beta3new(2);
f21 = beta3new(3);
f22 = beta3new(4);
f31 = beta3new(4);
f32 = beta3new(6);
f41 = beta3new(7);
f42 = beta3new(8);
f51 = beta3new(9);
f52 = beta3new(10);
f61 = beta3new(11);
f62 = beta3new(12);

n= length(dose);

nsize = wknum(:,2);

fvec1 = zeros([n 1]);
fvec1(1:sum(nsize(1))) = f11;
fvec1((sum(nsize(1))+1:sum(nsize(1:2)))) = f21;
fvec1((sum(nsize(1:2))+1:sum(nsize(1:3)))) = f31;
fvec1((sum(nsize(1:3))+1:sum(nsize(1:4)))) = f41;
fvec1((sum(nsize(1:4))+1:sum(nsize(1:5)))) = f51;
fvec1((sum(nsize(1:5))+1:sum(nsize(1:6)))) = f61;

fvec2 = zeros([n 1]);
fvec2(1:sum(nsize(1))) = f12;
fvec2((sum(nsize(1))+1:sum(nsize(1:2)))) = f22;
fvec2((sum(nsize(1:2))+1:sum(nsize(1:3)))) = f32;
fvec2((sum(nsize(1:3))+1:sum(nsize(1:4)))) = f42;
fvec2((sum(nsize(1:4))+1:sum(nsize(1:5)))) = f52;
fvec2((sum(nsize(1:5))+1:sum(nsize(1:6)))) = f62;
 


 difforignew = ((((fvec1)./( 1 + exp(slope1.*(dose - cen1)))) + ((1-fvec2-fvec1)./(1 + exp(slope2.*(dose - cen2)))) + ((fvec2)./(1 + exp(slope3.*(dose - cen3))))).*Vmaxall) - var;