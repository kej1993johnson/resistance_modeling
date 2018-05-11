function [new_vddata] = fit_by_var_by_time(vddata)
%This function will output a weight vector to be applied to weighting of
%all of the three different models. Instead of weighting based on the resi

% This script bins the dose axis into discrete sections that are determined
% by the data. They represent the large number of repeat doses at the same
% dose or very nearby.  For each bin, the average residual from the 1,2,
% and 3 population models are calculated and matched to the average of the
% dose value in that bin

% This will allow us to fit an equation to the average residual as a
% function of dose, which will then be used to weight the cost function so
% as to give higher weights to dose values with smaller residuals. The
% weight will be a function of dose.

% We will do this for each time point so that each week contains a
% different k


week = vddata(:,1);
dose = vddata (:,2);
var = vddata (:,3);
cohort_number = vddata(:,4);
n = length(dose);
time = length(unique(week));
vddata(:,5) = zeros([length(vddata) 1]);
week_numbers = unique(week);
for k = 1:time

    weekind = vddata(:,1) == week_numbers(k);
    oneweek = vddata(weekind,:); % creates vddata mx4 matrix that only contains data for one week


% 1st bin dose = 0
dosebin1 = [];
    for i=1:length(oneweek)
    if oneweek(i,2) < 0.0311
        dosebin1 = [dosebin1; oneweek(i,:)];
    end
    end
std_dosebin1 = std(dosebin1(:,3));
weight_to_add_1 = 1./std_dosebin1;
doseind1 = oneweek(:,2) < 0.0311;
oneweek(doseind1, 5) = weight_to_add_1;
 

 % 2nd bin dose = 0.03-1
dosebin2 = [];
    for i=1:length(oneweek)
    if oneweek(i,2) > 0.0311 && oneweek(i,2) <= 1
        dosebin2 = [dosebin2; oneweek(i,:)];
    end
    end
std_dosebin2 = std(dosebin2(:,3));
weight_to_add_2 = 1./std_dosebin2;
doseind2u = oneweek(:,2) <= 1;
doseind2l = oneweek(:,2) > 0.0311;
doseind2 = doseind2l==1 & doseind2u == 1;
oneweek(doseind2, 5) = weight_to_add_2;    
   
    
% 3rd bin dose = 1- 4
dosebin3 = [];
    for i=1:length(oneweek)
    if oneweek(i,2) > 1 && oneweek(i,2) <= 4
        dosebin3 = [dosebin3; oneweek(i,:)];
    end
    end
std_dosebin3 = std(dosebin3(:,3));
weight_to_add_3 = 1./std_dosebin3;
doseind3u = oneweek(:,2) <= 4;
doseind3l = oneweek(:,2) > 1;
doseind3 = doseind3l==1 & doseind3u == 1;
oneweek(doseind3, 5) = weight_to_add_3;

    

% 4th dose bin 4-8 w/8
dosebin4 = [];
    for i=1:length(oneweek)
    if oneweek(i,2) > 1 && oneweek(i,2) <= 4
        dosebin4 = [dosebin4; oneweek(i,:)];
    end
    end
std_dosebin4 = std(dosebin4(:,3));
weight_to_add_4 = 1./std_dosebin4;
doseind4u = oneweek(:,2) <= 8;
doseind4l = oneweek(:,2) > 4;
doseind4 = doseind4l==1 & doseind4u == 1;
oneweek(doseind4, 5) = weight_to_add_4;

% 5th dose bin 8-16 w/o 16
dosebin5 = [];
    for i=1:length(oneweek)
    if oneweek(i,2) > 8 && oneweek(i,2) < 16
        dosebin5 = [dosebin5; oneweek(i,:)];
    end
    end
std_dosebin5 = std(dosebin5(:,3));
weight_to_add_5 = 1./std_dosebin5;
doseind5u = oneweek(:,2) < 16;
doseind5l = oneweek(:,2) > 8;
doseind5 = doseind5l==1 & doseind5u == 1;
oneweek(doseind5, 5) = weight_to_add_5;

% 6th dose bin 16-24 w/16
dosebin6 = [];
    for i=1:length(oneweek)
    if oneweek(i,2) >= 16 && oneweek(i,2) < 24
        dosebin6 = [dosebin6; oneweek(i,:)];
    end
    end
std_dosebin6 = std(dosebin6(:,3));
weight_to_add_6 = 1./std_dosebin6;
doseind6u = oneweek(:,2) < 24;
doseind6l = oneweek(:,2) >= 16;
doseind6 = doseind6l==1 & doseind6u == 1;
oneweek(doseind6, 5) = weight_to_add_6;

 % 7th dose from 24-36 (not including 36)
 dosebin7 = [];
    for i=1:length(oneweek)
    if oneweek(i,2) >= 24 && oneweek(i,2) <= 36
        dosebin7 = [dosebin7; oneweek(i,:)];
    end
    end
std_dosebin7 = std(dosebin7(:,3));
weight_to_add_7 = 1./std_dosebin7;
doseind7u = oneweek(:,2) <= 36;
doseind7l = oneweek(:,2) >= 24;
doseind7 = doseind7l==1 & doseind7u == 1;
oneweek(doseind7, 5) = weight_to_add_7;
 

% 8th dose bin 36-45 
dosebin8 = [];
    for i=1:length(oneweek)
    if oneweek(i,2) > 36 && oneweek(i,2) < 45
        dosebin8 = [dosebin8; oneweek(i,:)];
    end
    end
std_dosebin8 = std(dosebin8(:,3));
weight_to_add_8 = 1./std_dosebin8;
doseind8u = oneweek(:,2) < 45;
doseind8l = oneweek(:,2) > 36;
doseind8 = doseind8l==1 & doseind8u == 1;
oneweek(doseind8, 5) = weight_to_add_8;

% 9th dose bine 45-54
dosebin9 = [];
    for i=1:length(oneweek)
    if oneweek(i,2) >= 45 && oneweek(i,2) < 54
        dosebin9 = [dosebin9; oneweek(i,:)];
    end
    end
if isempty(dosebin9)
else
std_dosebin9 = std(dosebin9(:,3));
weight_to_add_9 = 1./std_dosebin9;
doseind9u = oneweek(:,2) < 54;
doseind9l = oneweek(:,2) >= 45;
doseind9 = doseind9l==1 & doseind9u == 1;
oneweek(doseind9, 5) = weight_to_add_9;
end

% 10th dose bin 54-72
dosebin10 = [];
    for i=1:length(oneweek)
    if oneweek(i,2) >= 54 && oneweek(i,2) < 72
        dosebin10 = [dosebin10; oneweek(i,:)];
    end
    end
std_dosebin10 = std(dosebin10(:,3));
weight_to_add_10 = 1./std_dosebin10;
doseind10u = oneweek(:,2) < 72;
doseind10l = oneweek(:,2) >= 54;
doseind10 = doseind10l==1 & doseind10u == 1;
oneweek(doseind10, 5) = weight_to_add_10;

% 11th dose bin 72-96
dosebin11 = [];
    for i=1:length(oneweek)
    if oneweek(i,2) >= 72 && oneweek(i,2) < 96
        dosebin11 = [dosebin11; oneweek(i,:)];
    end
    end
std_dosebin11 = std(dosebin11(:,3));
weight_to_add_11 = 1./std_dosebin11;
doseind11u = oneweek(:,2) < 96;
doseind11l = oneweek(:,2) >= 72;
doseind11 = doseind11l==1 & doseind11u == 1;
oneweek(doseind11, 5) = weight_to_add_11; 

% 12th dose bin 96-109
dosebin12 = [];
    for i=1:length(oneweek)
    if oneweek(i,2) >= 96 && oneweek(i,2) < 109
        dosebin12 = [dosebin12; oneweek(i,:)];
    end
    end
if isempty(dosebin12)
else
std_dosebin12 = std(dosebin12(:,3));
weight_to_add_12 = 1./std_dosebin12;
doseind12u = oneweek(:,2) < 109;
doseind12l = oneweek(:,2) >= 96;
doseind12 = doseind12l==1 & doseind12u == 1;
oneweek(doseind12, 5) = weight_to_add_12;
end

% 13th bin 109-144
dosebin13 = [];
    for i=1:length(oneweek)
    if oneweek(i,2) >= 109 && oneweek(i,2) < 144
        dosebin13 = [dosebin13; oneweek(i,:)];
    end
    end
std_dosebin13 = std(dosebin13(:,3));
weight_to_add_13 = 1./std_dosebin13;
doseind13u = oneweek(:,2) < 144;
doseind13l = oneweek(:,2) >= 109;
doseind13 = doseind13l==1 & doseind13u == 1;
oneweek(doseind13, 5) = weight_to_add_13;

% 14th dose bin 144-146
dosebin14 = [];
    for i=1:length(oneweek)
    if oneweek(i,2) >= 144 && oneweek(i,2) < 146
        dosebin14 = [dosebin14; oneweek(i,:)];
    end
    end
std_dosebin14 = std(dosebin14(:,3));
weight_to_add_14 = 1./std_dosebin14;
doseind14u = oneweek(:,2) < 146;
doseind14l = oneweek(:,2) >= 144;
doseind14 = doseind14l==1 & doseind14u == 1;
oneweek(doseind14, 5) = weight_to_add_14;


% 15th dose bin 146-189
dosebin15 = [];
    for i=1:length(oneweek)
    if oneweek(i,2) >= 146 && oneweek(i,2) < 189
        dosebin15 = [dosebin15; oneweek(i,:)];
    end
    end
std_dosebin15 = std(dosebin15(:,3));
weight_to_add_15 = 1./std_dosebin15;
doseind15u = oneweek(:,2) < 189;
doseind15l = oneweek(:,2) >= 146;
doseind15 = doseind15l==1 & doseind15u == 1;
oneweek(doseind15, 5) = weight_to_add_15;


% 16th dose bin 189-220
dosebin16 = [];
    for i=1:length(oneweek)
    if oneweek(i,2) >= 189 && oneweek(i,2) < 220
        dosebin16 = [dosebin16; oneweek(i,:)];
    end
    end
std_dosebin16 = std(dosebin16(:,3));
weight_to_add_16 = 1./std_dosebin16;
doseind16u = oneweek(:,2) < 220;
doseind16l = oneweek(:,2) >= 189;
doseind16 = doseind16l==1 & doseind16u == 1;
oneweek(doseind16, 5) = weight_to_add_16;


onedoseind = oneweek(:,5) == Inf;
oneweek(onedoseind,5) = max(oneweek(isfinite(oneweek(:,5)),5));
cell_of_week_weights{k} = oneweek;

end
% Now just need to vertcat all the cells

new_vddata = [];
for k = 1:time
new_vddata = vertcat(new_vddata, cell2mat(cell_of_week_weights(k)));
end

largeind = new_vddata(:, 5) > 100;
    new_vddata(largeind,5) =50;
end

