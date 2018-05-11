function  v_model1allweeksnormed = model1popallweeksnormed( dose, betaLD50, Vmaxbyweek, nsize)
%This function produces a vector of of the data from the model


% Finds the average Vmax at each week and sets this as the numerator for
% the sigmoid function, changes as the code iterates through
num_weeks = length(Vmaxbyweek);
v_model1allweeksnormed= 1./( 1 + exp(betaLD50(1).*(dose(1:sum(nsize(1))) - betaLD50(num_weeks+1))));

 for i = 2: num_weeks
%v_model2 = Vmaxbyweek(i).*(beta2new(3+2*i)./( 1 + exp(beta2new(1).*(dose(sum(nsize(1:(i-1)))+1:sum(nsize(1:i))) - beta2new(2))))) + (beta2new(4+2*i)./(1 + exp(beta2new(3).*(dose(sum(nsize(1:(i-1)))+1:sum(nsize(1:i))) - beta2new(4)))));     
v_model1 = 1./( 1 + exp(betaLD50(i).*(dose(sum(nsize(1:(i-1)))+1:sum(nsize(1:i))) - betaLD50(num_weeks +i))));
v_model1allweeksnormed = vertcat(v_model1allweeksnormed, v_model1);
 end

end