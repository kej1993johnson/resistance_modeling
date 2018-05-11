function [weighted_diff_absresdens] = fit_simp1popabsresdens( beta1new, dose, var, k1, wknum, Vmaxall)
%This function is called by lsqnonlin.
% beta is a vector that contains the coefficeints of the equation. dose and
% var are the input data sets that were passed to lswqnonlin
slope1 = beta1new(1); % slope
cen1 = beta1new(2); % MD50 bulk

n= length(dose);

nsize = wknum(:,2);
 
difforignew = (Vmaxall./( 1 + exp(slope1.*(dose - cen1)))) - var;

weighted_diff_absres1= zeros([n,1]);
%weight_by_absres1 = %load('weight_by_absres1.mat');
weight_by_absres1 = wt_vector_resfit1( dose, k1 );
weight_by_absres1 = max(dose).*(weight_by_absres1./(sum(weight_by_absres1)));
%weight_by_absres1 = struct2cell(weight_by_absres1);
%weight_by_absres1 = cell2mat(weight_by_absres1);

for j = 1:n
weighted_diff_absres1(j) = weight_by_absres1(j)* difforignew(j);
end
weighted_diff_absres1 = weighted_diff_absres1';

% weighted_diff_density= zeros([n,1]);
% weight_by_density = load('weights_all_data_11_15.m');
% weight_by_density = load('weight_by_density.m')
%weight_by_density = (weight_by_density.*max(dose))./sum(weight_by_density);

% for k = 1:n
%     weighted_diff_density(k) = weight_by_density(k)* difforignew(k);
% end 
% weighted_diff_density = weighted_diff_density';


%weighted_diff_absresdens = 0.*weighted_diff_density + 1.*weighted_diff_absres1;

weighted_diff_absresdens = weighted_diff_absres1;

 
end