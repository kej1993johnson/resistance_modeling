function [ weights_by_fitted_res1] = wt_vector_resfit1( dose, k1 )
%This function creates a vector of weight values based on the model fitting
%of the absolute values of the residuals
model_fit1 = (k1(1)*(dose+k1(3)).*exp(-k1(2)*dose));
weights_by_fitted_res1 = 1./(model_fit1);

weights_by_fitted_res1 = (weights_by_fitted_res1 .* max(dose))./sum(weights_by_fitted_res1);

end

