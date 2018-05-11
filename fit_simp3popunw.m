function unweighted_diffnew = fit_simp3popunw( beta3new, dose, var, wknum)
%This function is called by lsqnonlin.
% beta is a vector that contains the coefficeints of the equation. dose and
% var are the input data sets that were passed to lswqnonlin
slope1 = beta3new(1); % slope 1 
cen1 = beta3new(2); % MD50 population 1
slope2 = beta3new(3);
cen2 = beta3new(4);
slope3 = beta3new(5);
cen3 =beta3new(6);
f11 = beta3new(7); % fraction of population 1 week 1
f12 = beta3new(8);
f13 = beta3new(9);
f21 = beta3new(10);
f22 = beta3new(11);
f23 = beta3new(12);
f31 = beta3new(13);
f32 = beta3new(14);
f33 = beta3new(15);
f41 = beta3new(16);
f42 = beta3new(17);
f43 = beta3new(18);
f51 = beta3new(19);
f52 = beta3new(20);
f53 = beta3new(21);
f61 = beta3new(22);
f62 = beta3new(23);
f63 = beta3new(24);
f71 = beta3new(25);
f72 = beta3new(26);
f73 = beta3new(27);
f81 = beta3new(28);
f82 = beta3new(29);
f83 = beta3new(30);
f91 = beta3new(31);
f92 = beta3new(32);
f93 = beta3new(33);
f101 = beta3new(34);
f102 = beta3new(35);
f103 = beta3new(36);
f111 = beta3new(37);
f112 = beta3new(38);
f113 = beta3new(39);
f121 = beta3new(40);
f122 = beta3new(41);
f123 = beta3new(42);


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

fvec2 = zeros( [n 1]);
fvec2(1:sum(nsize(1))) = f12;
fvec2((sum(nsize(1))+1:sum(nsize(1:2)))) = f22;
fvec2((sum(nsize(1:2))+1:sum(nsize(1:3)))) = f32;
fvec2((sum(nsize(1:3))+1:sum(nsize(1:4)))) = f42;
fvec2((sum(nsize(1:4))+1:sum(nsize(1:5)))) = f52;
fvec2((sum(nsize(1:5))+1:sum(nsize(1:6)))) = f62;
fvec2((sum(nsize(1:6))+1:sum(nsize(1:7)))) = f72;
fvec2((sum(nsize(1:7))+1:sum(nsize(1:8)))) = f82;
fvec2((sum(nsize(1:8))+1:sum(nsize(1:9)))) = f92;
fvec2((sum(nsize(1:9))+1:sum(nsize(1:10)))) = f102;
fvec2((sum(nsize(1:10))+1:sum(nsize(1:11)))) = f112;
fvec2((sum(nsize(1:11))+1:sum(nsize(1:12)))) = f122;

fvec3 = zeros( [n 1]);
fvec3(1:sum(nsize(1))) = f13;
fvec3((sum(nsize(1))+1:sum(nsize(1:2)))) = f23;
fvec3((sum(nsize(1:2))+1:sum(nsize(1:3)))) = f33;
fvec3((sum(nsize(1:3))+1:sum(nsize(1:4)))) = f43;
fvec3((sum(nsize(1:4))+1:sum(nsize(1:5)))) = f53;
fvec3((sum(nsize(1:5))+1:sum(nsize(1:6)))) = f63;
fvec3((sum(nsize(1:6))+1:sum(nsize(1:7)))) = f73;
fvec3((sum(nsize(1:7))+1:sum(nsize(1:8)))) = f83;
fvec3((sum(nsize(1:8))+1:sum(nsize(1:9)))) = f93;
fvec3((sum(nsize(1:9))+1:sum(nsize(1:10)))) = f103;
fvec3((sum(nsize(1:10))+1:sum(nsize(1:11)))) = f113;
fvec3((sum(nsize(1:11))+1:sum(nsize(1:12)))) = f123;

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




unweighted_diffnew = ((fvec1./(1 + exp(slope1.*(dose - cen1)))) + (fvec2./(1 + exp(slope2.*(dose - cen2)))) + (fvec3./(1 + exp(slope3*(dose - cen3)))).*Vmax) - var;


 
end