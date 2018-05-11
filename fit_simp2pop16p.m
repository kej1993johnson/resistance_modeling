function difforignew = fit_simp2pop16p( beta2new, dose, var )
%This function is called by lsqnonlin.
% beta is a vector that contains the coefficeints of the equation. dose and
% var are the input data sets that were passed to lswqnonlin
slope1 = beta2new(1); % slope
cen1 = beta2new(2); % MD50 bulk
slope2 = beta2new(3);
cen2 = beta2new(4);
f11 = beta2new(5); % fraction of population 1 week 1
f21 = beta2new(6); 
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

wknum = load ('week_size.m');
nsize = wknum(:,2);

n= length(dose);

%Find Vmax
Vmax1 = mean(var(1:12:nsize(1)));
Vmax2 = mean(var((nsize(1)+1):12:sum(nsize(1:2))));
Vmax3 = mean(var((sum(nsize(1:2))+1):12:sum(nsize(1:3))));
Vmax4 = mean(var((sum(nsize(1:3))+1):12:sum(nsize(1:4))));
Vmax5 = mean(var((sum(nsize(1:4))+1):12:sum(nsize(1:5))));
Vmax6 = mean(var((sum(nsize(1:5))+1):12:sum(nsize(1:6))));
Vmax7 = mean(var((sum(nsize(1:6))+1):12:sum(nsize(1:7))));
Vmax8 = mean(var((sum(nsize(1:7))+1):12:sum(nsize(1:8))));
Vmax9 = mean(var((sum(nsize(1:8))+1):12:sum(nsize(1:9))));
Vmax10 = mean(var((sum(nsize(1:9))+1):12:sum(nsize(1:10))));
Vmax11 = mean(var((sum(nsize(1:10))+1):12:sum(nsize(1:11))));
Vmax12 = mean(var((sum(nsize(1:11))+1):12:sum(nsize(1:12))));

%Replace Vmax
Vmax = zeros([n 1]);
Vmax(1:sum(nsize(1))) = Vmax1;
Vmax((sum(nsize(1))+1:sum(nsize(1:2)))) = Vmax2;
Vmax((sum(nsize(1:2))+1:sum(nsize(1:3)))) = Vmax3;
Vmax((sum(nsize(1:3))+1:sum(nsize(1:4)))) = Vmax4;
Vmax((sum(nsize(1:4))+1:sum(nsize(1:5)))) = Vmax5;
Vmax((sum(nsize(1:5))+1:sum(nsize(1:6)))) = Vmax6;
Vmax((sum(nsize(1:6))+1:sum(nsize(1:7)))) = Vmax7;
Vmax((sum(nsize(1:7))+1:sum(nsize(1:8)))) = Vmax8;
Vmax((sum(nsize(1:8))+1:sum(nsize(1:9)))) = Vmax9;
Vmax((sum(nsize(1:9))+1:sum(nsize(1:10)))) = Vmax10;
Vmax((sum(nsize(1:10))+1:sum(nsize(1:11)))) = Vmax11;
Vmax((sum(nsize(1:11))+1:sum(nsize(1:12)))) = Vmax12;



% totaldosefraction = zeros([n 1]);
%  for i = 2: n-1
%  
%  dosefraction = (dose(i+1) - dose(i-1))./(2*max(dose)); % tells you what fraction of the dose axis
%                                                         %this dose point
%                                                         %takes up
%  totaldosefraction (i,1) = dosefraction;
%   
% 
%  end 
%  totaldosefraction(1,1) = dose(2)./(max(dose));
%  totaldosefraction(n,1) =(dose(n) - dose(n-1))./(2*max(dose));
%  
wknum = load ('week_size.m');
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


 difforignew = ((Vmax).*(fvec1./( 1 + exp(slope1.*(dose - cen1)))) + ((1-fvec1)./(1 + exp(slope2.*(dose - cen2))))) - var;
 
% weighted_diffnew = zeros([n,1]);
% totaldosefraction = load('weights_all_data_10_5.m');
% for j = 1:n
% weighted_diffnew(j) = totaldosefraction(j)* difforignew(j);
% end
% weighted_diffnew = weighted_diffnew';

 
end