% This script loads the single treatment, repeat treatment at 2 and 6
% weeks, and untreated response. It combines fres and fsens parameters from
% the data fitting to estimate Nres and Nsens for each scenario



close all, clear all, clc

%% Run this code to get second treatment at 6 weeks analysis
% Load two population second treatment at 6 wks (16 wk total) data
% need to load proliferation data
fres2_twodose_6 = load('fres2twodoseall_6wks.mat');
error_barres2twodose_6 = load('errorbar2twodose_6wks.mat');
fres2_twodose_6 = struct2cell(fres2_twodose_6);
fres2_twodose_6 = cell2mat(fres2_twodose_6);
error_barres2twodose_6 = struct2cell(error_barres2twodose_6);
error_barres2twodose_6 = cell2mat(error_barres2twodose_6);

proliftwodose_6 = load('prolif_second_treatment_6.m');
proliftwodose_6reps = load('prolif_second_treatment_6_reps.m');

CItwodose_6_13 =  prctile(proliftwodose_6reps', [2.5 97.5]);
error_bar_two_dose6 = CItwodose_6_13(2,:) - CItwodose_6_13(1,:);

proliftwodose_6reps_2 = load('prolif_second_treatment_6_reps_2.m');

CItwodose_6_14_17 =  prctile(proliftwodose_6reps_2', [2.5 97.5]);
error_bar_two_dose6(1, 14:17) = CItwodose_6_14_17(2,:) - CItwodose_6_14_17(1,:);

t_wks = 15;
days = 7*(t_wks+1);
time = 0:1:(days-1); 

Ntot2_6 = zeros([days,1]);
Ntot2_6(1) = 1e6;
dt = 1;
fsens2_twodose_6 = 1-fres2_twodose_6;

Ntot_wks2_6(1) = 1e6;
Nres_wks2_6(1) = Ntot_wks2_6(1)*fres2_twodose_6(1);
for t = 2:t_wks+1
    Ntot_wks2_6(t) = Ntot_wks2_6(t-1).*proliftwodose_6(t).^7;
    Nres_wks2_6(t) = fres2_twodose_6(t)*Ntot_wks2_6(t);
    Nsens_wks2_6(t) = Ntot_wks2_6(t) - Nres_wks2_6(t);
end


%% Make into one column figure for second treatment at 6 weeks

subplot(3,1,2)
hold off
errorbar(0:1:16, proliftwodose_6, error_bar_two_dose6/2, '-ro', 'LineWidth',2);
xlabel ('Time Post Treatment (Weeks)', 'FontSize', 10)
ylabel ('Cells Per Day', 'FontSize', 10)
title('Twice Treated Proliferation', 'FontSize',10)
set(gca,'LineWidth',1.2,'FontSize',8)
ylim([0.5,1.6])
xlim([0 15])

subplot(3,1,1)
hold off
errorbar(0:1:15, fres2_twodose_6, error_barres2twodose_6./2, 'bo-', 'LineWidth', 1.2)
hold on
errorbar(0:1:15, fsens2_twodose_6, error_barres2twodose_6./2, 'yo-', 'LineWidth', 1.2)
xlim([0 15])
xlabel ('Time Post Treatment (Weeks)', 'FontSize', 10)
ylabel('Fraction of Cells','FontSize', 10)
title('Twice Treated Subpopulations','FontSize', 10)
set(gca,'LineWidth',1.2,'FontSize',8)
legend('Resistant', 'Sensitive')

subplot(3,1,3)
hold off
plot(0:1:15, log(Ntot_wks2_6),'go', 0:1:15, log(Nres_wks2_6), 'bo', 0:1:15, log(Nsens_wks2_6), 'yo', 'LineWidth', 2)
xlabel ('Time Post Treatment (Weeks)', 'FontSize', 10)
ylabel ('Log Cell Number)', 'FontSize', 10)
set(gca,'LineWidth',1.2,'FontSize',8)
title('Second Treatment at 6 weeks Phenotypic Dynamics','FontSize', 10)
legend('Total', 'Resistant', 'Sensitive')

%% Second Treatment at 2 Weeks
% Next load two population second treatment at 2 weeks (6 wks total) data
fres2_twodose_2 = load('fres2twodoseall_2wks.mat');
error_barres2twodose_2 = load('errorbar2twodose_2wks.mat');
fres2_twodose_2 = struct2cell(fres2_twodose_2);
fres2_twodose_2 = cell2mat(fres2_twodose_2);
error_barres2twodose_2 = struct2cell(error_barres2twodose_2);
error_barres2twodose_2 = cell2mat(error_barres2twodose_2);

proliftwodose_2 = load('prolif_second_treatment_2.m');

t_wks = 6;
days = 7*(t_wks+1);
time = 0:1:(days-1); 

dt = 1;
fsens2_twodose_2 = 1-fres2_twodose_2;

Ntot_wks2_2(1) = 1e6;
Nres_wks2_2(1) = Ntot_wks2_2(1)*fres2_twodose_2(1);
for t = 2:t_wks+1
    Ntot_wks2_2(t) = Ntot_wks2_2(t-1).*proliftwodose_2(t).^7;
    Nres_wks2_2(t) = fres2_twodose_2(t)*Ntot_wks2_2(t);
    Nsens_wks2_2(t) = Ntot_wks2_2(t) - Nres_wks2_2(t);
end



%% Make into one column figure for second treatment at 2 weeks

subplot(3,1,2)
hold off
plot(0:1:7, proliftwodose_2(1:8), '-ro', 'LineWidth',2);
xlabel ('Time Post Treatment (Weeks)', 'FontSize', 10)
ylabel ('Cells Per Day', 'FontSize', 10)
title('Twice Treated Proliferation', 'FontSize',10)
set(gca,'LineWidth',1.2,'FontSize',8)
ylim([0.5,1.6])
xlim([0 6])

subplot(3,1,1)
hold off
errorbar(0:1:6, fres2_twodose_2, error_barres2twodose_2./2, 'bo-', 'LineWidth', 1.2)
hold on
errorbar(0:1:6, fsens2_twodose_2, error_barres2twodose_2./2, 'yo-', 'LineWidth', 1.2)
xlim([0 6])
ylim([0 1])
xlabel ('Time Post Treatment (Weeks)', 'FontSize', 10)
ylabel('Fraction of Cells','FontSize', 10)
title('Twice Treated Subpopulations','FontSize', 10)
set(gca,'LineWidth',1.2,'FontSize',8)
legend('Resistant', 'Sensitive')

subplot(3,1,3)
hold off
plot(0:1:6, log(Ntot_wks2_2),'go', 0:1:6, log(Nres_wks2_2), 'bo', 0:1:6, log(Nsens_wks2_2), 'yo', 'LineWidth', 2)
xlabel ('Time Post Treatment (Weeks)', 'FontSize', 10)
ylabel ('Log Cell Number', 'FontSize', 10)
set(gca,'LineWidth',1.2,'FontSize',8)
title('Second Treatment at 2 weeks Phenotypic Dynamics','FontSize', 10)
legend('Total', 'Resistant', 'Sensitive')
%% Run this code to get single treatment 
load('fres2.mat');
prolif = load('prolif_single_treatment.m'); % in Number of cells per day
prolif = vertcat(prolif, prolif(end,1));

prolif_reps = load('prolif_single_treatment_reps.m');

CIsingledose_8 =  prctile(prolif_reps', [2.5 97.5]);
error_bar_single_dose = CIsingledose_8(2,:) - CIsingledose_8(1,:);

fsens2 = 1-fres2;

t_wks = 8;
time = 0:1: (days-1); 

Ntot(1) = 1e6;



Ntot_wks(1) = 1e6;
Nres_wks(1) = Ntot_wks(1)*fres2(1);
for t = 2:9
    Ntot_wks(t) = Ntot_wks(t-1).*prolif(t).^7;
    Nres_wks(t) = fres2(t)*Ntot_wks(t);
    Nsens_wks(t) = Ntot_wks(t) - Nres_wks(t);
end




errorbars_single_treat2_8 = load('errorbars_single_treat2_8.mat');
errorbars_single_treat2_8= struct2cell(errorbars_single_treat2_8);
errorbars_single_treat2_8 = cell2mat(errorbars_single_treat2_8);


%% Make into one column figure for single treatment
prolif_correct = prolif -1;
subplot(1,3,2)
hold off
errorbar(0:1:8, prolif_correct(1:9),error_bar_single_dose(1:9)/2, '-ro', 'LineWidth',2);
xlabel ('Time Post Treatment (Weeks)', 'FontSize', 10)
ylabel ('Drug Induced Proliferation Rate', 'FontSize', 10)
title('Single Treatment Proliferation', 'FontSize',10)
set(gca,'LineWidth',1.2,'FontSize',8)
ylim([-0.5,0.5])
xlim([0 8])

subplot(1,3,1)
hold off
errorbar(0:1:8, fres2, errorbars_single_treat2_8./2, 'bo-', 'LineWidth', 1.2)
hold on
errorbar(0:1:8, fsens2, errorbars_single_treat2_8./2, 'yo-', 'LineWidth', 1.2)
xlabel ('Time Post Treatment (Weeks)', 'FontSize', 10)
ylabel('Fraction of Cells','FontSize', 10)
title('Single Treatment Subpopulations','FontSize', 10)
set(gca,'LineWidth',1.2,'FontSize',8)
xlim([0 8])
legend('Resistant', 'Sensitive')
subplot(1,3,3)
hold off
plot(0:1:8, log(Ntot_wks),'g*', 0:1:8, log(Nres_wks), 'b*', 0:1:8, log(Nsens_wks), 'y*', 'LineWidth', 2)
xlabel ('Time Post Treatment (Weeks)', 'FontSize', 10)
ylabel ('Log Cell Number', 'FontSize', 10)
set(gca,'LineWidth',1.2,'FontSize',8)
xlim([0 8])
title('Single Treatment Phenotypic Dynamics','FontSize', 10)
legend('Total', 'Resistant', 'Sensitive')


%% Same analysis with untreated



load('fres2unt.mat');
prolif = load('prolif_untreated.m'); % in Number of cells per day
prolif = prolif(:,2);
guess_k = 1.33;
prolif = vertcat(prolif, guess_k, guess_k);
t_wks = 8;
days = 7*(t_wks+1);
time = 0:1: (days-1); 

Ntot_unt = zeros([days,1]);
Ntot_unt(1) = 1e6;
dt = 1;

prolif_tot_unt(1:7,1) = prolif(1).*ones([7,1]); % This will be ignored since it is proliferation before treatment

for i = 2:10
    
prolif_in_wk = prolif(i).*ones([7 1]);
prolif_tot_unt= vertcat(prolif_tot_unt, prolif_in_wk);
end

fres_unt_tot(1:7,1) = fres2unt(1).*ones([7,1]);
for i = 2:9   
fres_in_wk_unt = fres2unt(i).*ones([7 1]);
fres_unt_tot= vertcat(fres_unt_tot, fres_in_wk_unt);
end


for t = 2:days
Ntot_unt(t) = dt.*Ntot_unt(t-1).*prolif_tot_unt(t+7);
end

fsens2 = ones(days,1);
fsens_unt_tot = fsens2-fres_unt_tot;
Nres_unt = zeros(days,1);
Nsens_unt = zeros(days,1);
Nres_unt(1) = Ntot_unt(1)*fres_unt_tot(1);
for t=2:days
    Nres_unt(t)= fres_unt_tot(t)*Ntot_unt(t);
    Nsens_unt(t) = fsens_unt_tot(t)*Ntot_unt(t);
end


figure(2)

hold off
plot(1:1:70, prolif_tot_unt, '-o', 'LineWidth',2);
hold on
plot(1:1:70, prolif_correct, '-o', 'LineWidth',2);
xlabel ('Time Post Treatment (Days)', 'FontSize', 16)
ylabel ('Proliferation Rate (Number of Cells Per Day)', 'FontSize', 16)
title('Proliferation Rate vs. Time Post Treatment', 'FontSize',16)
legend('Untreated Daily Proliferation', 'Treated Daily Proliferation')
xlim([0,63])
ylim([0.5, 1.8])
set(gca,'LineWidth',1.5,'FontSize',16)

figure(3)
hold off
plot(time, fres_unt_tot, '-b', time, fsens_unt_tot,'-y', 'LineWidth', 2)
xlabel ('Time Post Treatment (Days)', 'FontSize', 16)
ylabel('Fraction of Cells','FontSize', 14)
title('Resistant and Sensitive Fraction Estimates','FontSize', 14)
set(gca,'LineWidth',1.5,'FontSize',16)
xlim([0,63])
legend('Resistant Fraction', 'Sensitive Fraction')
%%
subplot(2,2,1)
hold off
plot(time, fres_tot, '-b', time, fsens_tot,'-y', 'LineWidth', 2)
xlabel ('Time Post Treatment (Days)', 'FontSize', 14)
xlim([0 63])
ylim([0 1.2])
ylabel('Fraction of Cells','FontSize', 14)
title('Treated Population Fractions','FontSize', 14)
set(gca,'LineWidth',1.5,'FontSize',10)
legend('Resistant Fraction', 'Sensitive Fraction')


% %subplot(2,2,2)
% hold off
% plot(time, fres_unt_tot, '-b', time, fsens_unt_tot,'-y', 'LineWidth', 2)
% xlabel ('Time Post Treatment (Days)', 'FontSize', 14)
% xlim([0 63])
% ylim([0 1.2])
% ylabel('Fraction of Cells','FontSize', 14)
% title('Untreated Population Fractions','FontSize', 14)
% set(gca,'LineWidth',1.5,'FontSize',10)
% xlim([0,63])
% legend('Resistant Fraction', 'Sensitive Fraction')

subplot(3,2,2)
hold off
plot(1:1:70, prolif_tot_unt, '-o', 'LineWidth',2);
hold on
plot(1:1:70, prolif_tot, '-o', 'LineWidth',2);
xlabel ('Time Post Treatment (Days)', 'FontSize', 14)
ylabel ('Proliferation Rate (Number of Cells Per Day)', 'FontSize', 14)
title('Proliferation Rate vs. Time Post Treatment', 'FontSize',14)
legend('Untreated Daily Proliferation', 'Treated Daily Proliferation')
xlim([0,63])
ylim([0.5, 1.8])
set(gca,'LineWidth',1.5,'FontSize',10)

subplot(3,2,3)
hold off
plot(time, log(Ntot_unt),'g+', 'LineWidth', 2)
hold on
plot(time, log(Ntot),'go', 'LineWidth', 2)
xlabel ('Time Post Treatment (Days)', 'FontSize', 14)
ylabel ('Cell Number (log10 scale)', 'FontSize', 14)
set(gca,'LineWidth',1.5,'FontSize',8)
xlim([0 63])
title('Comparison of Treated and Untreated Cell Numbers','FontSize', 11)
legend('Total Untreated Cell Number ', 'Total Treated Cell Number')

subplot(2,2,4)
hold off
plot(time, log(Ntot),'go', time, log(Nres), 'bo', time, log(Nsens), 'yo', 'LineWidth', 2)
xlabel ('Time Post Treatment (Days)', 'FontSize', 14)
ylabel ('Cell Number (log10 scale)', 'FontSize', 14)
set(gca,'LineWidth',1.5,'FontSize',8)
xlim([0 63])
title('Model Estimate of Resistant and Sensitive Subpopulations','FontSize', 11)
legend('Total Cell Number', 'Number of Cells Resistant', 'Number of Cells Sensitive')

%% Make into three column figure

% first is single treatment
subplot(3,3,4)
hold off
errorbar(0:1:8, prolif(1:9),error_bar_single_dose(1:9)/2, '-ro', 'LineWidth',2);
xlabel ('Time Post Treatment (Weeks)', 'FontSize', 10)
ylabel ('Cells Per Day', 'FontSize', 10)
title('Single Treatment Proliferation', 'FontSize',10)
set(gca,'LineWidth',0.8,'FontSize',6)
ylim([0.5,1.6])
xlim([0,8])

subplot(3,3,1)
hold off
errorbar(0:1:8, fres2, errorbars_single_treat2_8./2, 'bo-', 'LineWidth', 1.2)
hold on
errorbar(0:1:8, fsens2, errorbars_single_treat2_8./2, 'yo-', 'LineWidth', 1.2)
xlabel ('Time Post Treatment (Weeks)', 'FontSize', 10)
ylabel('Fraction of Cells','FontSize', 10)
title('Single Treatment Subpopulations','FontSize', 8)
set(gca,'LineWidth',0.8,'FontSize',6)
legend('Resistant', 'Sensitive')
xlim([0,8])

subplot(3,3,7)
hold off
plot(0:1:8, log10(Ntot_wks),'g*', 0:1:8, log10(Nres_wks), 'b*', 0:1:8, log10(Nsens_wks), 'y*', 'LineWidth', 1.2)
xlabel ('Time Post Treatment (Weeks)', 'FontSize', 10)
ylabel ('Log Cell Number', 'FontSize', 10)
set(gca,'LineWidth',0.8,'FontSize',6)
title('Single Treatment Phenotypic Dynamics','FontSize', 8)
legend('Total', 'Resistant', 'Sensitive')
xlim([0,8])



%second is secodn treatment at 6 weeks
subplot(3,3,5)
hold off
errorbar(0:1:16, proliftwodose_6, error_bar_two_dose6/2, '-ro', 'LineWidth',2);
xlabel ('Time Post Treatment (Weeks)', 'FontSize', 10)
ylabel ('Cells Per Day', 'FontSize', 10)
title('Second Treatment at 6 Weeks Proliferation', 'FontSize',10)
set(gca,'LineWidth',0.8,'FontSize',6)
ylim([0.5,1.6])
xlim([0 15])

subplot(3,3,2)
hold off
errorbar(0:1:15, fres2_twodose_6, error_barres2twodose_6./2, 'bo-', 'LineWidth', 1.2)
hold on
errorbar(0:1:15, fsens2_twodose_6, error_barres2twodose_6./2, 'yo-', 'LineWidth', 1.2)
xlim([0 15])
xlabel ('Time Post Treatment (Weeks)', 'FontSize', 10)
ylabel('Fraction of Cells','FontSize', 10)
title('Second Treatment at 6 Weeks Subpopulations','FontSize', 8)
set(gca,'LineWidth',0.8,'FontSize',6)
%legend('Resistant', 'Sensitive')

subplot(3,3,8)
hold off
plot(0:1:15, log10(Ntot_wks2_6),'g*', 0:1:15, log10(Nres_wks2_6), 'b*', 0:1:15, log10(Nsens_wks2_6), 'y*', 'LineWidth', 2)
xlabel ('Time Post Treatment (Weeks)', 'FontSize', 10)
ylabel ('Log Cell Number', 'FontSize', 10)
set(gca,'LineWidth',0.8,'FontSize',6)
title('Second Treatment at 6 weeks Phenotypic Dynamics','FontSize', 8)
%legend('Total', 'Resistant', 'Sensitive')

% third is second treatment at 2 weeks

subplot(3,3,6)
hold off
plot(0:1:7, proliftwodose_2(1:8), '-ro', 'LineWidth',2);
xlabel ('Time Post Treatment (Weeks)', 'FontSize', 10)
ylabel ('Cells Per Day', 'FontSize', 10)
title('Second Treatment at 2 Weeks Proliferation', 'FontSize',10)
set(gca,'LineWidth',0.8,'FontSize',6)
ylim([0.5,1.6])
xlim([0 6])

subplot(3,3,3)
hold off
errorbar(0:1:6, fres2_twodose_2, error_barres2twodose_2./2, 'bo-', 'LineWidth', 1.2)
hold on
errorbar(0:1:6, fsens2_twodose_2, error_barres2twodose_2./2, 'yo-', 'LineWidth', 1.2)
xlim([0 6])
xlabel ('Time Post Treatment (Weeks)', 'FontSize', 10)
ylabel('Fraction of Cells','FontSize', 10)
title('Second Treatment at 2 Weeks Subpopulations','FontSize', 8)
set(gca,'LineWidth',0.8,'FontSize',6)
%legend('Resistant', 'Sensitive')
ylim([0 1])

subplot(3,3,9)
hold off
plot(0:1:6, log10(Ntot_wks2_2),'g*', 0:1:6, log10(Nres_wks2_2), 'b*', 0:1:6, log10(Nsens_wks2_2), 'y*', 'LineWidth', 1.2)
xlabel ('Time Post Treatment (Weeks)', 'FontSize', 10)
ylabel ('Log Cell Number', 'FontSize', 10)
set(gca,'LineWidth',0.8,'FontSize',6)
title('Second Treatment at 2 Weeks Phenotypic Dynamics','FontSize', 8)
%legend('Total', 'Resistant', 'Sensitive')
