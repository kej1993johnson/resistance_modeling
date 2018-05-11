function [ v_model1 ] = model1pop( dose, beta1new, Vmaxavg )
%This function produces a vector of of the data from the model
%Find Vmax



v_model1 = (1./( 1 + exp(beta1new(1).*(dose - beta1new(2))))) ;


end

