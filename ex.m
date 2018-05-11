
model_bulk_LD50 = Vmaxbyweek(1)./(1 + exp(betaLD50_8(1).*(X - betaLD50_8((1+(length(betaLD50_8)./2))))));
%%
for i = 2: length(Vmaxbyweek)
model_bulk_LD50_new2 = Vmaxbyweek(i)./(1 + exp(betaLD50_8(i).*(X - betaLD50_8((i+(length(betaLD50_8)./2)))))); 
model_bulk_LD50 = vertcat(model_bulk_LD50, model_bulk_LD50_new2);
end