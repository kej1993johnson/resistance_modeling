function [ v_model2 ] = model2pop12p( dose, beta2new, wk, Vmaxbyweek)
%This function produces a vector of of the data from the model
Vmax = repmat(Vmaxbyweek(wk), length(dose),1);

v_model2 = Vmax.*((beta2new(4+wk)./( 1 + exp(beta2new(1).*(dose - beta2new(2))))) + (1-beta2new(4+wk))./(1 + exp(beta2new(3).*(dose - beta2new(4)))));
% = Vmax.*(beta2new(4+wk)./( 1 + exp(beta2new(1).*(dose - beta2new(2))))) + ((1-beta2new(4+wk))./(1 + exp(beta2new(3).*(dose - beta2new(4))))); 


end