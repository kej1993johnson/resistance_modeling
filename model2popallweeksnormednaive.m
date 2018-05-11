function  v_model2allweeksnormednaive = model2popallweeksnormednaive( dosenc, coeffs2naive, beta2naive_8)
%This function produces a vector of of the data from the model


% Finds the average Vmax at each week and sets this as the numerator for
% the sigmoid function, changes as the code iterates through

v_model2allweeksnormednaive = (beta2naive_8./( 1 + exp(coeffs2naive(1).*(dosenc - coeffs2naive(2))))) + ((1-beta2naive_8)./(1 + exp(coeffs2naive(3).*(dosenc - coeffs2naive(4)))));



end