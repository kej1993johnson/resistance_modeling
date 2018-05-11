function  v_model3allweeksnormed = model3popallweeksnormed( dose, beta3new, Vmaxbyweek, nsize)
%This function produces a vector of of the data from the model


% Finds the average Vmax at each week and sets this as the numerator for
% the sigmoid function, changes as the code iterates through

v_model3allweeksnormed= (beta3new(7)./( 1 + exp(beta3new(1).*(dose(1:sum(nsize(1))) - beta3new(2))))) + ((beta3new(8))./(1 + exp(beta3new(3).*(dose(1:sum(nsize(1))) - beta3new(4))))) + ((1-beta3new(7) - beta3new(8))./(1 + exp(beta3new(5).*(dose(1:sum(nsize(1))) - beta3new(6)))));


 for i = 2: length(Vmaxbyweek)
%v_model2 = Vmaxbyweek(i).*(beta2new(3+2*i)./( 1 + exp(beta2new(1).*(dose(sum(nsize(1:(i-1)))+1:sum(nsize(1:i))) - beta2new(2))))) + (beta2new(4+2*i)./(1 + exp(beta2new(3).*(dose(sum(nsize(1:(i-1)))+1:sum(nsize(1:i))) - beta2new(4)))));     
v_model3 = (beta3new(5+2*i)./( 1 + exp(beta3new(1).*(dose(sum(nsize(1:(i-1)))+1:sum(nsize(1:i))) - beta3new(2)))))+(beta3new(6+2*i)./( 1 + exp(beta3new(3).*(dose(sum(nsize(1:(i-1)))+1:sum(nsize(1:i))) - beta3new(4))))) + ((1-beta3new(6+2*i)-beta3new(5+2*i))./(1 + exp(beta3new(5).*(dose(sum(nsize(1:(i-1)))+1:sum(nsize(1:i))) - beta3new(8)))));
v_model3allweeksnormed = vertcat(v_model3allweeksnormed, v_model3);
 end

end