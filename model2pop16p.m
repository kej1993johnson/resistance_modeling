function [ v_model2 ] = model2pop16p( dose, beta2new, Vmaxbyweek, wk )
%This function produces a vector of of the data from the model
v_model2 = Vmaxbyweek(wk,1).*(beta2new(4+wk)./( 1 + exp(beta2new(1).*(dose - beta2new(2))))) + (1-beta2new(4+wk))./(1 + exp(beta2new(3).*(dose - beta2new(4))));
%v_model2 = Vmax.*(beta2new(4+wk)./( 1 + exp(beta2new(1).*(dose - beta2new(2))))) + ((1-beta2new(4+wk))./(1 + exp(beta2new(3).*(dose - beta2new(4))))); 


end