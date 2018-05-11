function  v_model2allweeksnormed = model2popallweeksnormed( dose, beta2new, Vmaxbyweek, nsize)
%This function produces a vector of of the data from the model


% Finds the average Vmax at each week and sets this as the numerator for
% the sigmoid function, changes as the code iterates through

v_model2allweeksnormed= 1.*(beta2new(5)./( 1 + exp(beta2new(1).*(dose(1:sum(nsize(1))) - beta2new(2))))) + ((1-beta2new(5))./(1 + exp(beta2new(3).*(dose(1:sum(nsize(1))) - beta2new(4)))));

 for i = 2: length(Vmaxbyweek)
%v_model2 = Vmaxbyweek(i).*(beta2new(3+2*i)./( 1 + exp(beta2new(1).*(dose(sum(nsize(1:(i-1)))+1:sum(nsize(1:i))) - beta2new(2))))) + (beta2new(4+2*i)./(1 + exp(beta2new(3).*(dose(sum(nsize(1:(i-1)))+1:sum(nsize(1:i))) - beta2new(4)))));     
v_model2 = 1.*(beta2new(4+i)./( 1 + exp(beta2new(1).*(dose(sum(nsize(1:(i-1)))+1:sum(nsize(1:i))) - beta2new(2))))) + (1-beta2new(4+i))./(1 + exp(beta2new(3).*(dose(sum(nsize(1:(i-1)))+1:sum(nsize(1:i))) - beta2new(4))));
v_model2allweeksnormed = vertcat(v_model2allweeksnormed, v_model2);
 end

end