function difforignew = fit_simp3pop6weeksnaive( beta3new, dose, var, coeffs3 )
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



 


 difforignew = ((((f11)./( 1 + exp(slope1.*(dose - cen1)))) + ((1-f12-f11)./(1 + exp(slope2.*(dose - cen2)))) + ((f12)./(1 + exp(slope3.*(dose - cen3))))).*var(1)) - var;