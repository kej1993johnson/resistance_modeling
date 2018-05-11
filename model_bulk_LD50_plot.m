function [ model_bulk_LD50 ] = model_bulk_LD50_plot( X, betaLD50, Vmaxbyweek)


%This function produces a vector of of the data from the model

model_bulk_LD50 = Vmaxbyweek(1)./(1 + exp(betaLD50(1).*(X - betaLD50((1+(length(betaLD50)./2))))));

for i = 2: length(Vmaxbyweek)
model_bulk_LD50_new2 = Vmaxbyweek(i)./(1 + exp(betaLD50(i).*(X - betaLD50((i+(length(betaLD50)./2)))))); 
model_bulk_LD50 = vertcat(model_bulk_LD50, model_bulk_LD50_new2);
end

end