function  v_model3allweeksnormednaive = model3popallweeksnormednaive( dose, coeffs3, beta3naive)
%This function produces a vector of of the data from the model


% Finds the average Vmax at each week and sets this as the numerator for
% the sigmoid function, changes as the code iterates through

v_model3allweeksnormednaive= (beta3naive(1)./( 1 + exp(coeffs3(1).*(dose - coeffs3(2))))) + ((beta3naive(2))./(1 + exp(coeffs3(3).*(dose - coeffs3(4))))) + ((1-beta3naive(1) - beta3naive(2))./(1 + exp(coeffs3(5).*(dose - coeffs3(6)))));




end