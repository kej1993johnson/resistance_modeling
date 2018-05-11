function [ model_bulk_LD50 ] = model_bulk_LD50( dose, betaLD50, Vmaxall, wk )

% Rewrite this function so that it changes by week but so that model is all
% weeks
%This function produces a vector of of the data from the model
model_bulk_LD50 = Vmaxall./(1 + exp(betaLD50(wk).*(dose - betaLD50((wk+(length(betaLD50)./2)))))); 
end
