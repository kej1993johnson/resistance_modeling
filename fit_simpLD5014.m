function difforignew = fit_simpLD5014( betaLD50, dose, var, wknum, Vmaxall )
%This function is called by lsqnonlin. It fits all of the data to a
%different slope and LD50 value at each week.  

slope1 = betaLD50(1); 
slope2 = betaLD50(2); 
slope3 = betaLD50(3);
slope4 = betaLD50(4);
slope5 = betaLD50(5);
slope6 = betaLD50(6);
slope7 = betaLD50(7);
slope8 = betaLD50(8);
slope9 = betaLD50(9);
slope10 = betaLD50(10);
slope11 = betaLD50(11);
slope12 = betaLD50(12);
slope13 = betaLD50(13);
slope14 = betaLD50(14);
cen1 = betaLD50(15);
cen2 = betaLD50(16);
cen3 = betaLD50(17);
cen4 = betaLD50(18);
cen5 = betaLD50(19);
cen6 = betaLD50(20);
cen7 = betaLD50(21);
cen8 = betaLD50(22);
cen9 = betaLD50(23);
cen10 = betaLD50(24);
cen11 = betaLD50(25);
cen12 = betaLD50(26);
cen13 = betaLD50(27);
cen14 = betaLD50(28);



n= length(dose);

nsize = wknum(:,2);

slope = zeros([n 1]);
slope(1:sum(nsize(1))) = slope1;
slope((sum(nsize(1))+1:sum(nsize(1:2)))) = slope2;
slope((sum(nsize(1:2))+1:sum(nsize(1:3)))) = slope3;
slope((sum(nsize(1:3))+1:sum(nsize(1:4)))) = slope4;
slope((sum(nsize(1:4))+1:sum(nsize(1:5)))) = slope5;
slope((sum(nsize(1:5))+1:sum(nsize(1:6)))) = slope6;
slope((sum(nsize(1:6))+1:sum(nsize(1:7)))) = slope7;
slope((sum(nsize(1:7))+1:sum(nsize(1:8)))) = slope8;
slope((sum(nsize(1:8))+1:sum(nsize(1:9)))) = slope9;
slope((sum(nsize(1:9))+1:sum(nsize(1:10)))) = slope10;
slope((sum(nsize(1:10))+1:sum(nsize(1:11)))) = slope11;
slope((sum(nsize(1:11))+1:sum(nsize(1:12)))) = slope12;
slope((sum(nsize(1:12))+1:sum(nsize(1:13)))) = slope13;
slope((sum(nsize(1:13))+1:sum(nsize(1:14)))) = slope14;

cen = zeros([n 1]);
cen(1:sum(nsize(1))) = cen1;
cen((sum(nsize(1))+1:sum(nsize(1:2)))) = cen2;
cen((sum(nsize(1:2))+1:sum(nsize(1:3)))) = cen3;
cen((sum(nsize(1:3))+1:sum(nsize(1:4)))) = cen4;
cen((sum(nsize(1:4))+1:sum(nsize(1:5)))) = cen5;
cen((sum(nsize(1:5))+1:sum(nsize(1:6)))) = cen6;
cen((sum(nsize(1:6))+1:sum(nsize(1:7)))) = cen7;
cen((sum(nsize(1:7))+1:sum(nsize(1:8)))) = cen8;
cen((sum(nsize(1:8))+1:sum(nsize(1:9)))) = cen9;
cen((sum(nsize(1:9))+1:sum(nsize(1:10)))) = cen10;
cen((sum(nsize(1:10))+1:sum(nsize(1:11)))) = cen11;
cen((sum(nsize(1:11))+1:sum(nsize(1:12)))) = cen12;
cen((sum(nsize(1:12))+1:sum(nsize(1:13)))) = cen13;
cen((sum(nsize(1:13))+1:sum(nsize(1:14)))) = cen14;



 difforignew = Vmaxall./( 1 + exp(slope.*(dose - cen))) - var;