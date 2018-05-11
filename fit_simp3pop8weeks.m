function difforignew = fit_simp3pop8weeks( beta3new, dose, var, wknum, Vmaxall )
%This function is called by lsqnonlin.
% beta is a vector that contains the coefficeints of the equation. dose and
% var are the input data sets that were passed to lswqnonlin
slope1 = beta3new(1); % slope
cen1 = beta3new(2); % MD50 bulk
slope2 = beta3new(3);
cen2 = beta3new(4);
slope3 = beta3new(5);
cen3 = beta3new(6);
f11 = beta3new(7);
f12 = beta3new(8);
f21 = beta3new(9);
f22 = beta3new(10);
f31 = beta3new(11);
f32 = beta3new(12);
f41 = beta3new(13);
f42 = beta3new(14);
f51 = beta3new(15);
f52 = beta3new(16);
f61 = beta3new(17);
f62 = beta3new(18);
f71 = beta3new(19);
f72 = beta3new(20);
f81 = beta3new(21);
f82 = beta3new(22);

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

fvec2 = zeros([n 1]);
fvec2(1:sum(nsize(1))) = f12;
fvec2((sum(nsize(1))+1:sum(nsize(1:2)))) = f22;
fvec2((sum(nsize(1:2))+1:sum(nsize(1:3)))) = f32;
fvec2((sum(nsize(1:3))+1:sum(nsize(1:4)))) = f42;
fvec2((sum(nsize(1:4))+1:sum(nsize(1:5)))) = f52;
fvec2((sum(nsize(1:5))+1:sum(nsize(1:6)))) = f62;
fvec2((sum(nsize(1:6))+1:sum(nsize(1:7)))) = f72;
fvec2((sum(nsize(1:7))+1:sum(nsize(1:8)))) = f82;


 difforignew = ((((fvec1)./( 1 + exp(slope1.*(dose - cen1)))) + ((fvec2)./(1 + exp(slope2.*(dose - cen2)))) + ((1-fvec1-fvec2)./(1 + exp(slope3.*(dose - cen3))))).*Vmaxall) - var;