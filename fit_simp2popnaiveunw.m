function difforignew = fit_simp2popnaiveunw( beta2new, dose, var, coeffs2, Vmaxnaiveavg )
%This function is called by lsqnonlin.beta2new, dose, var, k2,wknum, Vmaxall )
% beta is a vector that contains the coefficeints of the equation. dose and
% var are the input data sets that were passed to lswqnonlin
f11 = beta2new(1); % fraction of population 1 week 0


n = length(dose);



difforignew = Vmaxnaiveavg.*(var(1).*(((f11))./( 1 + exp(coeffs2(1).*(dose - coeffs2(2))))) + ((1-f11)./(1 + exp(coeffs2(3).*(dose - coeffs2(4)))))) - var;


 
end