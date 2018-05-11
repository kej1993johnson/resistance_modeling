function [k1,k2] = fitabsoluteresiduals( dose, var, wknum, Vmaxall)
% This script bins the dose axis into discrete sections that are determined
% by the data. They represent the large number of repeat doses at the same
% dose or very nearby.  For each bin, the average residual from the 1,2,
% and 3 population models are calculated and matched to the average of the
% dose value in that bin

% This will allow us to fit an equation to the average residual as a
% function of dose, which will then be used to weight the cost function so
% as to give higher weights to dose values with smaller residuals. The
% weight will be a function of dose.
 
n = length(dose);

% FINDS unweighted residual vectors based on unweighted model function
initials1new = [ .05; 30]; 
[beta1newunw, resnorm1unw, residuals1unw] = lsqnonlin(@fit_simp1popunw, initials1new,[0; 0],[ Inf; Inf],[],dose,var, wknum, Vmaxall);

options = optimset('Display','off','FunValCheck','on', ...
                   'MaxFunEvals',Inf,'MaxIter',Inf, ...
                   'TolFun',1e-6,'TolX',1e-6);

paramslb = zeros([16 1]);
paramsub = [ Inf; Inf; Inf; Inf; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1];
params0 = [ .1; 17; .04; 50; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5; .5];
[beta2newunw, resnorm2unw, residuals2unw] = lsqnonlin(@fit_simp2popunw16p,...
    params0,...
    paramslb,...
    paramsub,...
    options,...
    dose,...
    var,...
    wknum,...
    Vmaxall);

initials3 = [ .2; 13; .07; 30; .05; 60; .25; .4; .35; .25; .4; .35; .25; .4; .35; .25; .4; .35; .25; .4; .35; .25; .4; .35; .25; .4; .35; .25; .4; .35; .25; .4; .35; .25; .4; .35; .25; .4; .35; .25; .4; .35];
paramslb3 = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;];
paramsub3 = [ Inf; Inf; Inf; Inf; Inf; Inf; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1];

[beta3newunw, resnorm3unw, residuals3unw] = lsqnonlin(@fit_simp3popunw,...
    initials3,...
    paramslb3,...
    paramsub3,...
    options,...
    dose,...
    var,...
    wknum);
% Next take their absolute value
absresiduals1 = abs(residuals1unw);
absresiduals2 = abs(residuals2unw);
absresiduals3 = abs(residuals3unw);

% create a matrix of n rows with columns (dose, res1, res2, res3)
dose_res = horzcat(dose,absresiduals1, absresiduals2, absresiduals3);

% creates matrix with 3 rows of zeros to be filled by weights for the
% specific dose
data_w_weights = [dose_res zeros([ n 3])]; 

% 1st bin dose = 0
% This section iterates through the data and finds all 0 doses (so doses
% less than 0.0311 uM). For that row, it extracts the dose, and the
% residual values of the model at that row.  Then it puts these new arrays
% corresponding to the first bin, finds the averages of each value, and
% then puts the averages back into a 1X4 matrix containing avg dose and avg
% abs(residual)for each model


bin1 = {};
bin1res1 = {};
bin1res2 = {};
bin1res3 = {};
for i=1:n
    if dose_res(i,1) < 0.0311
        bin1 = [bin1; dose_res(i,1)];
        bin1res1 = [bin1res1; dose_res(i,2)];
        bin1res2 = [bin1res2; dose_res(i,3)];
        bin1res3 = [bin1res3; dose_res(i,4)];
    end
end
bin1 = cell2mat(bin1);
bin1res1 = cell2mat(bin1res1);
bin1res2 = cell2mat(bin1res2);
bin1res3 = cell2mat(bin1res3);
ndose1 = length(bin1);
avgdose1= mean(bin1);
avgbin1res1 = mean(bin1res1);
avgbin1res2 = mean(bin1res2);
avgbin1res3 = mean(bin1res3);
% returns the average dose and residual for each model
% now want to put into dose, avg_res matrix

dose_resavg1 = horzcat(avgdose1, avgbin1res1, avgbin1res2, avgbin1res3);
% Later, will take this and say dos_resavgallbins = vertcat(dose_resavg1,
% dos_resavg2, etc...)

w1_1 = 1./(avgbin1res1);
w1_2 = 1./(avgbin1res2);
w1_3 = 1./(avgbin1res3);

for j = 1:n
if dose_res(j,1) < 0.0311
    data_w_weights(j,5) = w1_1;
    data_w_weights(j,6) = w1_2;
    data_w_weights(j,7) = w1_3;
end
end
% puts the weight for each row of the data into the 5, 6, or 7 columns


% 2nd bin dose = 0.03-1

bin2 = {};
bin2res1 = {};
bin2res2 = {};
bin2res3 = {};
for i=1:n
    if dose_res(i,1)> 0.0311 && dose_res(i,1) <= 1
        bin2 = [bin2; dose_res(i,1)];
        bin2res1 = [bin2res1; dose_res(i,2)];
        bin2res2 = [bin2res2; dose_res(i,3)];
        bin2res3 = [bin2res3; dose_res(i,4)];
    end
end
bin2 = cell2mat(bin2);
bin2res1 = cell2mat(bin2res1);
bin2res2 = cell2mat(bin2res2);
bin2res3 = cell2mat(bin2res3);
ndose2 = length(bin2);
avgdose2= mean(bin2);
avgbin2res1 = mean(bin2res1);
avgbin2res2 = mean(bin2res2);
avgbin2res3 = mean(bin2res3);

% returns the average dose and residual for each model
% now want to put into dose, avg_res matrix

dose_resavg2 = horzcat(avgdose2, avgbin2res1, avgbin2res2, avgbin2res3);

w2_1 = 1./(avgbin2res1);
w2_2 = 1./(avgbin2res2);
w2_3 = 1./(avgbin2res3);

for j = 1:n
if dose_res(j,1)> 0.0311 && dose_res(j,1) <= 1
    data_w_weights(j,5) = w2_1;
    data_w_weights(j,6) = w2_2;
    data_w_weights(j,7) = w2_3;
end
end

% 3rd bin dose = 1- 4

bin3 = {};
bin3res1 = {};
bin3res2 = {};
bin3res3 = {};
for i=1:n
    if dose_res(i,1) > 1 && dose_res(i,1) <= 4 
        bin3 = [bin3; dose_res(i,1)];
        bin3res1 = [bin3res1; dose_res(i,2)];
        bin3res2 = [bin3res2; dose_res(i,3)];
        bin3res3 = [bin3res3; dose_res(i,4)];
    end
end
bin3 = cell2mat(bin3);
bin3res1 = cell2mat(bin3res1);
bin3res2 = cell2mat(bin3res2);
bin3res3 = cell2mat(bin3res3);
ndose3 = length(bin3);
avgdose3= mean(bin3);
avgbin3res1 = mean(bin3res1);
avgbin3res2 = mean(bin3res2);
avgbin3res3 = mean(bin3res3);

% returns the average dose and residual for each model
% now want to put into dose, avg_res matrix

dose_resavg3 = horzcat(avgdose3, avgbin3res1, avgbin3res2, avgbin3res3);

w3_1 = 1./(avgbin3res1);
w3_2 = 1./(avgbin3res2);
w3_3 = 1./(avgbin3res3);

for j = 1:n
    if dose_res(j,1) > 1 && dose_res(j,1) <= 4  
    data_w_weights(j,5) = w3_1;
    data_w_weights(j,6) = w3_2;
    data_w_weights(j,7) = w3_3;
    end
end

% 4th dose bin 4-8 w/8
bin4 = {};
bin4res1 = {};
bin4res2 = {};
bin4res3 = {};
for i=1:n
    if dose_res(i,1) >= 4 && dose_res(i,1) <= 8 
        bin4 = [bin4; dose_res(i,1)];
        bin4res1 = [bin4res1; dose_res(i,2)];
        bin4res2 = [bin4res2; dose_res(i,3)];
        bin4res3 = [bin4res3; dose_res(i,4)];
    end
end
bin4 = cell2mat(bin4);
bin4res1 = cell2mat(bin4res1);
bin4res2 = cell2mat(bin4res2);
bin4res3 = cell2mat(bin4res3);
ndose4 = length(bin4);
avgdose4= mean(bin4);
avgbin4res1 = mean(bin4res1);
avgbin4res2 = mean(bin4res2);
avgbin4res3 = mean(bin4res3);

% returns the average dose and residual for each model
% now want to put into dose, avg_res matrix

dose_resavg4 = horzcat(avgdose4, avgbin4res1, avgbin4res2, avgbin4res3);

w4_1 = 1./(avgbin4res1);
w4_2 = 1./(avgbin4res2);
w4_3 = 1./(avgbin4res3);

for j = 1:n
    if dose_res(j,1) >= 4 && dose_res(j,1) <= 8 
    data_w_weights(j,5) = w4_1;
    data_w_weights(j,6) = w4_2;
    data_w_weights(j,7) = w4_3;
    end
end

% 5th dose bin 8-16 w/o 16
bin5 = {};
bin5res1 = {};
bin5res2 = {};
bin5res3 = {};
for i=1:n
    if dose_res(i,1) > 8 && dose_res(i,1) < 16
        bin5 = [bin5; dose_res(i,1)];
        bin5res1 = [bin5res1; dose_res(i,2)];
        bin5res2 = [bin5res2; dose_res(i,3)];
        bin5res3 = [bin5res3; dose_res(i,4)];
    end
end
bin5 = cell2mat(bin5);
bin5res1 = cell2mat(bin5res1);
bin5res2 = cell2mat(bin5res2);
bin5res3 = cell2mat(bin5res3);
ndose5 = length(bin5);
avgdose5= mean(bin5);
avgbin5res1 = mean(bin5res1);
avgbin5res2 = mean(bin5res2);
avgbin5res3 = mean(bin5res3);

% returns the average dose and residual for each model
% now want to put into dose, avg_res matrix

dose_resavg5 = horzcat(avgdose5, avgbin5res1, avgbin5res2, avgbin5res3);

w5_1 = 1./(avgbin5res1);
w5_2 = 1./(avgbin5res2);
w5_3 = 1./(avgbin5res3);

for j = 1:n
    if dose_res(j,1) > 8 && dose_res(j,1) < 16 
    data_w_weights(j,5) = w5_1;
    data_w_weights(j,6) = w5_2;
    data_w_weights(j,7) = w5_3;
    end
end

% 6th dose bin 16-24 w/16

bin6 = {};
bin6res1 = {};
bin6res2 = {};
bin6res3 = {};
for i=1:n
    if dose_res(i,1) >=16 && dose_res(i,1) < 24
        bin6 = [bin6; dose_res(i,1)];
        bin6res1 = [bin6res1; dose_res(i,2)];
        bin6res2 = [bin6res2; dose_res(i,3)];
        bin6res3 = [bin6res3; dose_res(i,4)];
    end
end
bin6 = cell2mat(bin6);
bin6res1 = cell2mat(bin6res1);
bin6res2 = cell2mat(bin6res2);
bin6res3 = cell2mat(bin6res3);
ndose6 = length(bin6);
avgdose6= mean(bin6);
avgbin6res1 = mean(bin5res1);
avgbin6res2 = mean(bin6res2);
avgbin6res3 = mean(bin6res3);

% returns the average dose and residual for each model
% now want to put into dose, avg_res matrix

dose_resavg6 = horzcat(avgdose6, avgbin6res1, avgbin6res2, avgbin6res3);

w6_1 = 1./(avgbin6res1);
w6_2 = 1./(avgbin6res2);
w6_3 = 1./(avgbin6res3);

for j = 1:n
    if dose_res(j,1) >=16 && dose_res(j,1) < 24 
    data_w_weights(j,5) = w6_1;
    data_w_weights(j,6) = w6_2;
    data_w_weights(j,7) = w6_3;
    end
end

 % 7th dose from 24-36 (not including 36)
bin7 = {};
bin7res1 = {};
bin7res2 = {};
bin7res3 = {};
for i=1:n
    if dose_res(i,1) >=24 && dose_res(i,1) < 36
        bin7 = [bin7; dose_res(i,1)];
        bin7res1 = [bin7res1; dose_res(i,2)];
        bin7res2 = [bin7res2; dose_res(i,3)];
        bin7res3 = [bin7res3; dose_res(i,4)];
    end
end
bin7 = cell2mat(bin7);
bin7res1 = cell2mat(bin7res1);
bin7res2 = cell2mat(bin7res2);
bin7res3 = cell2mat(bin7res3);
ndose7 = length(bin7);
avgdose7= mean(bin7);
avgbin7res1 = mean(bin7res1);
avgbin7res2 = mean(bin7res2);
avgbin7res3 = mean(bin7res3);

% returns the average dose and residual for each model
% now want to put into dose, avg_res matrix

dose_resavg7 = horzcat(avgdose7, avgbin7res1, avgbin7res2, avgbin7res3);

w7_1 = 1./(avgbin7res1);
w7_2 = 1./(avgbin7res2);
w7_3 = 1./(avgbin7res3);

for j = 1:n
    if dose_res(j,1) >=24 && dose_res(j,1) < 36 
    data_w_weights(j,5) = w7_1;
    data_w_weights(j,6) = w7_2;
    data_w_weights(j,7) = w7_3;
    end
end

% 8th dose bin 36-45 
bin8 = {};
bin8res1 = {};
bin8res2 = {};
bin8res3 = {};
for i=1:n
    if dose_res(i,1) >=36 && dose_res(i,1) < 45
        bin8 = [bin8; dose_res(i,1)];
        bin8res1 = [bin8res1; dose_res(i,2)];
        bin8res2 = [bin8res2; dose_res(i,3)];
        bin8res3 = [bin8res3; dose_res(i,4)];
    end
end
bin8 = cell2mat(bin8);
bin8res1 = cell2mat(bin8res1);
bin8res2 = cell2mat(bin8res2);
bin8res3 = cell2mat(bin8res3);
ndose8 = length(bin8);
avgdose8= mean(bin8);
avgbin8res1 = mean(bin8res1);
avgbin8res2 = mean(bin8res2);
avgbin8res3 = mean(bin8res3);

% returns the average dose and residual for each model
% now want to put into dose, avg_res matrix

dose_resavg8 = horzcat(avgdose8, avgbin8res1, avgbin8res2, avgbin8res3);
w8_1 = 1./(avgbin8res1);
w8_2 = 1./(avgbin8res2);
w8_3 = 1./(avgbin8res3);

for j = 1:n
    if dose_res(j,1) >=36 && dose_res(j,1) < 45
    data_w_weights(j,5) = w8_1;
    data_w_weights(j,6) = w8_2;
    data_w_weights(j,7) = w8_3;
    end
end

% 9th dose bine 45-54
bin9 = {};
bin9res1 = {};
bin9res2 = {};
bin9res3 = {};
for i=1:n
    if dose_res(i,1) >=45 && dose_res(i,1) < 54
        bin9 = [bin9; dose_res(i,1)];
        bin9res1 = [bin9res1; dose_res(i,2)];
        bin9res2 = [bin9res2; dose_res(i,3)];
        bin9res3 = [bin9res3; dose_res(i,4)];
    end
end
bin9 = cell2mat(bin9);
bin9res1 = cell2mat(bin9res1);
bin9res2 = cell2mat(bin9res2);
bin9res3 = cell2mat(bin9res3);
ndose9 = length(bin9);
avgdose9= mean(bin9);
avgbin9res1 = mean(bin9res1);
avgbin9res2 = mean(bin9res2);
avgbin9res3 = mean(bin9res3);

% returns the average dose and residual for each model
% now want to put into dose, avg_res matrix

dose_resavg9 = horzcat(avgdose9, avgbin9res1, avgbin9res2, avgbin9res3);

w9_1 = 1./(avgbin9res1);
w9_2 = 1./(avgbin9res2);
w9_3 = 1./(avgbin9res3);

for j = 1:n
    if dose_res(j,1) >=45 && dose_res(j,1) < 54
    data_w_weights(j,5) = w9_1;
    data_w_weights(j,6) = w9_2;
    data_w_weights(j,7) = w9_3;
    end
end


% 10th dose bin 54-72

bin10 = {};
bin10res1 = {};
bin10res2 = {};
bin10res3 = {};
for i=1:n
    if dose_res(i,1) >=54 && dose_res(i,1) < 72
        bin10 = [bin10; dose_res(i,1)];
        bin10res1 = [bin10res1; dose_res(i,2)];
        bin10res2 = [bin10res2; dose_res(i,3)];
        bin10res3 = [bin10res3; dose_res(i,4)];
    end
end
bin10 = cell2mat(bin10);
bin10res1 = cell2mat(bin10res1);
bin10res2 = cell2mat(bin10res2);
bin10res3 = cell2mat(bin10res3);
ndose10 = length(bin10);
avgdose10= mean(bin10);
avgbin10res1 = mean(bin10res1);
avgbin10res2 = mean(bin10res2);
avgbin10res3 = mean(bin10res3);

% returns the average dose and residual for each model
% now want to put into dose, avg_res matrix

dose_resavg10 = horzcat(avgdose10, avgbin10res1, avgbin10res2, avgbin10res3);

w10_1 = 1./(avgbin10res1);
w10_2 = 1./(avgbin10res2);
w10_3 = 1./(avgbin10res3);

for j = 1:n
    if dose_res(j,1) >=54 && dose_res(j,1) < 72
    data_w_weights(j,5) = w10_1;
    data_w_weights(j,6) = w10_2;
    data_w_weights(j,7) = w10_3;
    end
end

% 11th dose bin 72-96
 
bin11 = {};
bin11res1 = {};
bin11res2 = {};
bin11res3 = {};
for i=1:n
    if dose_res(i,1) >=72 && dose_res(i,1) < 96
        bin11 = [bin11; dose_res(i,1)];
        bin11res1 = [bin11res1; dose_res(i,2)];
        bin11res2 = [bin11res2; dose_res(i,3)];
        bin11res3 = [bin11res3; dose_res(i,4)];
    end
end
bin11 = cell2mat(bin11);
bin11res1 = cell2mat(bin11res1);
bin11res2 = cell2mat(bin11res2);
bin11res3 = cell2mat(bin11res3);
ndose11 = length(bin11);
avgdose11= mean(bin11);
avgbin11res1 = mean(bin11res1);
avgbin11res2 = mean(bin11res2);
avgbin11res3 = mean(bin11res3);

% returns the average dose and residual for each model
% now want to put into dose, avg_res matrix

dose_resavg11 = horzcat(avgdose11, avgbin11res1, avgbin11res2, avgbin11res3);

w11_1 = 1./(avgbin11res1);
w11_2 = 1./(avgbin11res2);
w11_3 = 1./(avgbin11res3);

for j = 1:n
    if dose_res(j,1) >=72 && dose_res(j,1) < 96
    data_w_weights(j,5) = w11_1;
    data_w_weights(j,6) = w11_2;
    data_w_weights(j,7) = w11_3;
    end
end

% 12th dose bin 96-109

bin12 = {};
bin12res1 = {};
bin12res2 = {};
bin12res3 = {};
for i=1:n
    if dose_res(i,1) >=96 && dose_res(i,1) < 109
        bin12 = [bin12; dose_res(i,1)];
        bin12res1 = [bin12res1; dose_res(i,2)];
        bin12res2 = [bin12res2; dose_res(i,3)];
        bin12res3 = [bin12res3; dose_res(i,4)];
    end
end
bin12 = cell2mat(bin12);
bin12res1 = cell2mat(bin12res1);
bin12res2 = cell2mat(bin12res2);
bin12res3 = cell2mat(bin12res3);
ndose12 = length(bin12);
avgdose12= mean(bin12);
avgbin12res1 = mean(bin12res1);
avgbin12res2 = mean(bin12res2);
avgbin12res3 = mean(bin12res3);

% returns the average dose and residual for each model
% now want to put into dose, avg_res matrix

dose_resavg12 = horzcat(avgdose12, avgbin12res1, avgbin12res2, avgbin12res3);

w12_1 = 1./(avgbin12res1);
w12_2 = 1./(avgbin12res2);
w12_3 = 1./(avgbin12res3);

for j = 1:n
    if dose_res(j,1) >=96 && dose_res(j,1) < 109
    data_w_weights(j,5) = w12_1;
    data_w_weights(j,6) = w12_2;
    data_w_weights(j,7) = w12_3;
    end
end
% 13th bin 109-144

bin13 = {};
bin13res1 = {};
bin13res2 = {};
bin13res3 = {};
for i=1:n
    if dose_res(i,1) >=109 && dose_res(i,1) < 144
       bin13 = [bin13; dose_res(i,1)];
        bin13res1 = [bin13res1; dose_res(i,2)];
        bin13res2 = [bin13res2; dose_res(i,3)];
        bin13res3 = [bin13res3; dose_res(i,4)];
    end
end
bin13 = cell2mat(bin13);
bin13res1 = cell2mat(bin13res1);
bin13res2 = cell2mat(bin13res2);
bin13res3 = cell2mat(bin13res3);
ndose13 = length(bin13);
avgdose13= mean(bin13);
avgbin13res1 = mean(bin13res1);
avgbin13res2 = mean(bin13res2);
avgbin13res3 = mean(bin13res3);

% returns the average dose and residual for each model
% now want to put into dose, avg_res matrix

dose_resavg13 = horzcat(avgdose13, avgbin13res1, avgbin13res2, avgbin13res3);

w13_1 = 1./(avgbin13res1);
w13_2 = 1./(avgbin13res2);
w13_3 = 1./(avgbin13res3);

for j = 1:n
    if dose_res(j,1) >=109 && dose_res(j,1) < 144
    data_w_weights(j,5) = w13_1;
    data_w_weights(j,6) = w13_2;
    data_w_weights(j,7) = w13_3;
    end
end

% 14th dose bin 144-146
bin14 = {};
bin14res1 = {};
bin14res2 = {};
bin14res3 = {};
for i=1:n
    if dose_res(i,1) >=144 && dose_res(i,1) < 146
        bin14 = [bin14; dose_res(i,1)];
        bin14res1 = [bin14res1; dose_res(i,2)];
        bin14res2 = [bin14res2; dose_res(i,3)];
        bin14res3 = [bin14res3; dose_res(i,4)];
    end
end
bin14 = cell2mat(bin14);
bin14res1 = cell2mat(bin14res1);
bin14res2 = cell2mat(bin14res2);
bin14res3 = cell2mat(bin14res3);
ndose14 = length(bin14);
avgdose14= mean(bin14);
avgbin14res1 = mean(bin14res1);
avgbin14res2 = mean(bin14res2);
avgbin14res3 = mean(bin14res3);

% returns the average dose and residual for each model
% now want to put into dose, avg_res matrix

dose_resavg14 = horzcat(avgdose14, avgbin14res1, avgbin14res2, avgbin14res3);

w14_1 = 1./(avgbin14res1);
w14_2 = 1./(avgbin14res2);
w14_3 = 1./(avgbin14res3);

for j = 1:n
    if dose_res(j,1) >=144 && dose_res(j,1) < 146
    data_w_weights(j,5) = w14_1;
    data_w_weights(j,6) = w14_2;
    data_w_weights(j,7) = w14_3;
    end
end

% 15th dose bin 146-189
bin15 = {};
bin15res1 = {};
bin15res2 = {};
bin15res3 = {};
for i=1:n
    if dose_res(i,1) >=146 && dose_res(i,1) < 189
        bin15 = [bin15; dose_res(i,1)];
        bin15res1 = [bin15res1; dose_res(i,2)];
        bin15res2 = [bin15res2; dose_res(i,3)];
        bin15res3 = [bin15res3; dose_res(i,4)];
    end
end
bin15 = cell2mat(bin15);
bin15res1 = cell2mat(bin15res1);
bin15res2 = cell2mat(bin15res2);
bin15res3 = cell2mat(bin15res3);
ndose15 = length(bin15);
avgdose15= mean(bin15);
avgbin15res1 = mean(bin15res1);
avgbin15res2 = mean(bin15res2);
avgbin15res3 = mean(bin15res3);

% returns the average dose and residual for each model
% now want to put into dose, avg_res matrix

dose_resavg15 = horzcat(avgdose15, avgbin15res1, avgbin15res2, avgbin15res3);

w15_1 = 1./(avgbin15res1);
w15_2 = 1./(avgbin15res2);
w15_3 = 1./(avgbin15res3);

for j = 1:n
    if dose_res(j,1) >=146 && dose_res(j,1) < 189
    data_w_weights(j,5) = w15_1;
    data_w_weights(j,6) = w15_2;
    data_w_weights(j,7) = w15_3;
    end
end

% 16th dose bin 189-220
bin16 = {};
bin16res1 = {};
bin16res2 = {};
bin16res3 = {};
for i=1:n
    if dose_res(i,1) >=189 && dose_res(i,1) < 220
        bin16 = [bin16; dose_res(i,1)];
        bin16res1 = [bin16res1; dose_res(i,2)];
        bin16res2 = [bin16res2; dose_res(i,3)];
        bin16res3 = [bin16res3; dose_res(i,4)];
    end
end
bin16 = cell2mat(bin16);
bin16res1 = cell2mat(bin16res1);
bin16res2 = cell2mat(bin16res2);
bin16res3 = cell2mat(bin16res3);
ndose16 = length(bin16);
avgdose16= mean(bin16);
avgbin16res1 = mean(bin16res1);
avgbin16res2 = mean(bin16res2);
avgbin16res3 = mean(bin16res3);
% returns the average dose and residual for each model
% now want to put into dose, avg_res matrix

dose_resavg16 = horzcat(avgdose16, avgbin16res1, avgbin16res2, avgbin16res2);

w16_1 = 1./(avgbin16res1);
w16_2 = 1./(avgbin16res2);
w16_3 = 1./(avgbin16res3);

for j = 1:n
    if dose_res(j,1) >=189 && dose_res(j,1) < 220
    data_w_weights(j,5) = w16_1;
    data_w_weights(j,6) = w16_2;
    data_w_weights(j,7) = w16_3;
    end
end
%Create one matrix of avg dose and avg residual at that dose
dose_res_avg_all = vertcat(dose_resavg1, dose_resavg2, dose_resavg3, dose_resavg4,dose_resavg5, dose_resavg6, dose_resavg7, dose_resavg8,dose_resavg9, dose_resavg10, dose_resavg11, dose_resavg12,dose_resavg13, dose_resavg14, dose_resavg15, dose_resavg16);


%
% figure (1)
% hold off
% plot(dose_res_avg_all(:,1), dose_res_avg_all(:,2), 'o')
% xlabel('dose (uM')
% ylabel('absolute value of residual')
% title(' Average Abolsute Value of Residual vs. Dose One Population Model')
% 
% 
% figure (2)
% hold off
% plot(dose_res_avg_all(:,1), dose_res_avg_all(:,3), 'o')
% xlabel('dose (uM')
% ylabel('absolute value of residual')
% title(' Average Abolsute Value of Residual vs. Dose Two Population Model')
% 
% figure (3)
% hold off
% plot(dose_res_avg_all(:,1), dose_res_avg_all(:,4), 'o')
% xlabel('dose (uM')
% ylabel('absolute value of residual')
% title(' Average Abolsute Value of Residual vs. Dose Three Population Model')

% Now Create and Save Weight Vector

weight_by_absres1 = data_w_weights(:,5);
weight_by_absres2 = data_w_weights(:,6);
weight_by_absres3 = data_w_weights(:,7);



% Fit Absolute values of residuals of bulk model

doseavg = dose_res_avg_all(:,1);
res1avg = dose_res_avg_all(:,2);

k1guess = [0.01; .03; 15];

[k1, newresnorm1, newres1] = lsqnonlin(@fitabsres1, k1guess,[0; 0; 0],[ Inf; Inf; Inf],[],doseavg,res1avg);

k1 = k1;
%exp(slope1.*(dose - cen1)))
X = (.1:.1:250);
model_fit1 = (k1(1)*(X+k1(3)).*exp(-k1(2)*X));
W1 = 1./(model_fit1);

% figure(1)
% hold off
% plot( doseavg, res1avg, 'o');
% hold on
% plot(X, model_fit1,'LineWidth',4);
% xlim ([0, 250])
% xlabel('dose (uM)','FontSize',16)
% ylabel('absolute value of residual','FontSize',16)
% title('Fitted Function for Absolute Value of Residual for One Population Model','FontSize',15)
% 
% figure(2)
% hold off
% plot(X, W1, 'LineWidth', 4)
% xlim([0,250])
% xlabel('dose(uM)', 'FontSize',16)
% ylabel('weight', 'FontSize',16)
% title('Weighting Scheme for One Population Model', 'FontSize',15)


% Fit absolute values of residuals of two population model

res2avg = dose_res_avg_all(:,3);

k2guess = [0.01; .03; 7];

[k2, newresnorm2, newres2] = lsqnonlin(@fitabsres2, k2guess,[0; 0; 0],[ Inf; Inf; Inf],[],doseavg,res2avg);

k2 = k2;

X = (0:.1:250);
model_fit2 = (k2(1)*(X+k2(3)).*exp(-k2(2)*X));
W2 = 1./(model_fit2);

% figure(3)
% hold off
% plot( doseavg, res2avg, 'o');
% hold on
% plot(X, model_fit2, 'LineWidth', 4);
% xlabel('dose (uM)','FontSize',16)
% ylabel('absolute value of residual','FontSize',16)
% title('Fitted Function for Absolute Value of Residual Two Population Model','FontSize',15)
% 
% figure(4)
% hold off
% plot(X, W2, 'LineWidth', 4)
% xlim([0,250])
% xlabel('dose(uM)', 'FontSize',16)
% ylabel('weight', 'FontSize',16)
% title('Weighting Scheme for Two Population Model', 'FontSize',15)
% Fit absolute values of residuals of three population model

res3avg = dose_res_avg_all(:,4);

k3guess = [0.01; .03; 7];

[k3, newresnorm3, newres3] = lsqnonlin(@fitabsres3, k3guess,[0; 0; 0],[ Inf; Inf; Inf],[],doseavg,res3avg);

k3 = k3;

X = (0:.1:250);
model_fit3 = (k3(1)*(X+k3(3)).*exp(-k3(2)*X));
W3 = 1./model_fit3;

% figure(5)
% hold off
% plot( doseavg, res3avg, 'o');
% hold on
% plot(X, model_fit3, 'LineWidth', 4);
% xlabel('dose (uM)','FontSize',16)
% ylabel('absolute value of residual','FontSize',16)
% title('Fitted Function for Absolute Value of Residual Three Population Model','FontSize',15)
% 
% figure (6)
% hold off
% plot(X, W3, 'LineWidth', 4)
% xlim([0,250])
% xlabel('dose(uM)', 'FontSize',16)
% ylabel('weight', 'FontSize',16)
% title('Weighting Scheme for Three Population Model', 'FontSize',15)


end

