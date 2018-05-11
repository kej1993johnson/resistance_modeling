function weighted_diffabsresdens = fit_simp2popabsresdensnormed( beta2new, dose, var, k2,wknum, Vmaxall )
%This function is called by lsqnonlin.
% beta is a vector that contains the coefficeints of the equation. dose and
% var are the input data sets that were passed to lswqnonlin
slope1 = beta2new(1); % slope
cen1 = beta2new(2); % MD50 bulk
slope2 = beta2new(3);
cen2 = beta2new(4);
f11 = beta2new(5);
f21 = beta2new(6);% fraction of population 1 week 1
f31 = beta2new(7);
f41 = beta2new(8);
f51 = beta2new(9);
f61 = beta2new(10);
f71 = beta2new(11);
f81 = beta2new(12);
f91 = beta2new(13);
f101 = beta2new(14);
f111 = beta2new(15);
f121 = beta2new(16);


n= length(dose);

nsize = wknum(:,2);

fvec1 = zeros([n 1]);
fvec1(1:sum(nsize(1))) = f11;
fvec1((sum(nsize(1))+1:sum(nsize(1:2)))) = f21;
fvec1((sum(nsize(1:2))+1:sum(nsize(1:3)))) = f31;
fvec1((sum(nsize(1:3))+1:sum(nsize(1:4)))) = f41;
fvec1((sum(nsize(1:4))+1:sum(nsize(1:5)))) = f51;
fvec1((sum(nsize(1:5))+1:sum(nsize(1:6)))) = f61;
fvec1((sum(nsize(1:6))+1:sum(nsize(1:7)))) = f71;
fvec1((sum(nsize(1:7))+1:sum(nsize(1:8)))) = f81;
fvec1((sum(nsize(1:8))+1:sum(nsize(1:9)))) = f91;
fvec1((sum(nsize(1:9))+1:sum(nsize(1:10)))) = f101;
fvec1((sum(nsize(1:10))+1:sum(nsize(1:11)))) = f111;
fvec1((sum(nsize(1:11))+1:sum(nsize(1:12)))) = f121;
    
 difforignew = ((((fvec1)./( 1 + exp(slope1.*(dose - cen1)))) + ((1-fvec1)./(1 + exp(slope2.*(dose - cen2)))) + 0).*Vmaxall) - var;
 
weighted_diffabsres2= zeros([n,1]);
weight_by_absres2 = wt_vector_resfit1( dose, k2 );
weight_by_absres2 = max(dose).*(weight_by_absres2./(sum(weight_by_absres2)));

for j = 1:n
weighted_diffabsres2(j) = weight_by_absres2(j)* difforignew(j);
end
weighted_diffabsres2 = weighted_diffabsres2';

% weighted_diff_density= zeros([n,1]);
% weight_by_density = load('weights_all_data_11_15.m');
% weight_by_density = (weight_by_density.*max(dose))./sum(weight_by_density);
% 
% 
% for k = 1:n
%     weighted_diff_density(k) = weight_by_density(k)* difforignew(k);
% end 
% weighted_diff_density = weighted_diff_density';


% weighted_diffabsresdens = 0.*weighted_diff_density + 1.*weighted_diffabsres2;
weighted_diffabsresdens = weighted_diffabsres2;

 
end