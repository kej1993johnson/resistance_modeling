function difforignew = fit_simp2popunw16p( beta2new, dose, var, wknum, Vmaxall )
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

% fvec2 = zeros( [n 1]);
% fvec2(1:sum(nsize(1))) = f12;
% fvec2((sum(nsize(1))+1:sum(nsize(1:2)))) = f22;
% fvec2((sum(nsize(1:2))+1:sum(nsize(1:3)))) = f32;
% fvec2((sum(nsize(1:3))+1:sum(nsize(1:4)))) = f42;
% fvec2((sum(nsize(1:4))+1:sum(nsize(1:5)))) = f52;
% fvec2((sum(nsize(1:5))+1:sum(nsize(1:6)))) = f62;
% fvec2((sum(nsize(1:6))+1:sum(nsize(1:7)))) = f72;
% fvec2((sum(nsize(1:7))+1:sum(nsize(1:8)))) = f82;
% fvec2((sum(nsize(1:8))+1:sum(nsize(1:9)))) = f92;
% fvec2((sum(nsize(1:9))+1:sum(nsize(1:10)))) = f102;
% fvec2((sum(nsize(1:10))+1:sum(nsize(1:11)))) = f112;
% fvec2((sum(nsize(1:11))+1:sum(nsize(1:12)))) = f122;

% Finds the average Vmax at each week and sets this as the numerator for
% the sigmoid function, changes as the code iterates through

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




    

 difforignew = ((((1-fvec1)./( 1 + exp(slope1.*(dose - cen1)))) + ((fvec1)./(1 + exp(slope2.*(dose - cen2)))) + 0).*Vmaxall) - var;