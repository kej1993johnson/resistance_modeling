

% Figures cleaned up (ish)

figure (1);
hold off
subplot(1, 3, 1)
% single static population model
hold off
plot (dose, var, 'g.', 'LineWidth', 1.5)
hold on
plot(dose,var_unt, 'k.', 'LineWidth', 1.5)
plot (X, v_model1plot, '-g', 'LineWidth', 3);
plot(X, v_model1plotunt, '-k', 'LineWidth',3);
xlim([0,220]);
xlabel('Dose (\muM)','FontSize',10)
ylabel('Viability','FontSize',10)
%title ('Single Static Population Model','FontSize',18)
legend('Treated LD50 = 50.4 +/-2.4 \muM', 'Untreated LD50 = 37.0 +/- 3.5 \muM')
legend boxoff
set(gca,'LineWidth',1.5,'FontSize',10)
title('a                        ', 'FontSize', 18)

subplot(1, 3, 2)
hold off
% single dynamic population model
errorbar(time, LD50s, errorbarsLD50./2, '-go', 'LineWidth', 1.5)
hold on
xlabel ('Time (weeks)', 'FontSize', 10)
ylabel ('LD50 (\muM)', 'FontSize', 10)
set(gca,'LineWidth',1.5,'FontSize',10)
errorbar(time, LD50sunt, errorbarsLD50unt./2,'ko', 'LineWidth',1.5)
ylim ([20,90])
xlim([-0.2 8.2])
title('b                        ', 'FontSize', 18)
legend('Treated LD50', 'Untreated LD50')
legend boxoff

subplot(1,3, 3)
% two population model
hold off
errorbar(time, fres2, errorbars./2, 'ro-', 'LineWidth', 1.5)
hold on
errorbar(time, fsens2, errorbars./2, 'bo-', 'LineWidth', 1.5)
ylim([ 0 1])
xlim([-0.2 8.2])
set(gca,'LineWidth',1.5,'FontSize',10)
xlabel('Time (weeks)', 'FontSize', 10)
ylabel('Fraction of Cells', 'FontSize',10)
%title ('Two Population Model', 'FontSize', 14)
legend ('Resistant LD50 = 79.7 +/- 6.5  \muM', 'Sensitive LD50 = 22.4 +/- 2.0 \muM')
legend boxoff
title('c                        ', 'FontSize', 18)
%%
figure (2)

subplot(1,2,1)
% example of weeks 2 and 8 fit
hold off
plot (dose(sum(nsize(1:1))+1:sum(nsize(1:2))), var(sum(nsize(1:1))+1:sum(nsize(1:2))), 's', 'LineWidth', 0.6)
hold on
plot (dose(sum(nsize(1:7))+1:sum(nsize(1:8))), var(sum(nsize(1:7))+1:sum(nsize(1:8))), '>', 'LineWidth', 0.6)
plot (X, v_model2_2, '-r', 'LineWidth', 3);
plot (X, v_model2_8, '-b', 'LineWidth', 3);
legend( 'Data 2 weeks', 'Data 8 weeks', 'Model 2 weeks', 'Model 8 weeks')
xlim ([0 220])
xlabel('dose (\muM)','FontSize',11)
ylabel('Viability','FontSize',11)
title ('a','FontSize',18)
legend boxoff

fracs = [1-beta2newall(6), beta2newall(6); 1-beta2newall(12), beta2newall(12)];
errorbarsfracs = [errorbars(3), errorbars(3); errorbars(9), errorbars(9)];
ctrs = 1:2;
data = fracs;

c = categorical({'8 weeks_{res}', '8 weeks_{sens}'});
subplot(1,2,2)
hold off
hb = bar(fracs, 'grouped')
hold on
% For each set of bars, find the centers of the bars, and write error bars
pause(0.1); %pause allows the figure to be created
for ib = 1:numel(hb)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hb(ib).XData+hb(ib).XOffset;
    errorbar(xData,fracs(:,ib),0.5.*errorbarsfracs(:,ib),'k.')
end
hold on
%errorbar( fracs, errorbarsfracs)
title ('b','FontSize',18)
ylabel ('Fraction of Cells')
legend ('resistant', 'sensitive')
legend boxoff
% bar graph of weeks 2 and 8

%%
figure(7)
hBar = bar(ctrs,data);
hold on


for k1= 1:length(hBar)
    hb = get(hBar(k1),'XData');
    midbar = mean(hb);
    errorbar(midbar, data(:,k1), errorbarsfracs(:,k1), '.')
    sigbarx(k1,:) = midbar;
end
%     ctr(k1,:) = bsxfun(@plus,hBar(1).Xdata, [hBar(k1).Xoffset]');
%     ydt(k1,:) = hBar(k1).Ydata;
% end
% hold on
% errorbar(ctr', ydt', errorbarsfracs, '.r')
% hold off

%
%% Code gets proliferation data and two population resistant and sensitive cell data 
 % figure 3


load('fres2.mat');
% prolif = load('prolif_single_treatment.m'); % in Number of cells per day
% prolif = vertcat(prolif, prolif(end,1));

prolif_reps = load('prolif_single_treatment_reps.m');
prolif_mean = mean(prolif_reps, 2) % gets mean at each time point (1st is pretreatment, last is at 9 WPT
%
CIsingledose_8 =  prctile(prolif_reps', [2.5 97.5]);
error_bar_single_dose = CIsingledose_8(2,:) - CIsingledose_8(1,:);

fsens2 = 1-fres2;

t_wks = 8;

Ntot(1) = 1e6;
errorbars_single_treat2_8 = load('errorbars_single_treat2_8.mat');
errorbars_single_treat2_8= struct2cell(errorbars_single_treat2_8);
errorbars_single_treat2_8 = cell2mat(errorbars_single_treat2_8);

freslow = fres2- errorbars_single_treat2_8./2;
freshigh = fres2 + errorbars_single_treat2_8./2;

Ntot_wks(1) = 1e6;
Nres_wks(1) = Ntot_wks(1)*fres2(1);
Nsens_wks(1) = Ntot_wks(1)-Nres_wks(1);
for t = 2:9
    Ntot_wks(t) = Ntot_wks(t-1).*prolif_mean(t).^7;
    Nres_wks(t) = fres2(t)*Ntot_wks(t);
    Nsens_wks(t) = Ntot_wks(t) - Nres_wks(t);
end

% Find upper and lower CI

Nlowtot(1) = 1e6;
Nhightot(1) = 1e6;

for t = 2:9
    Nlowtot(t) = Ntot_wks(t-1).*CIsingledose_8(1,t).^7;
    Nhightot(t) = Ntot_wks(t-1).*CIsingledose_8(2,t).^7;
    Nreslow(t)= freslow(t)*Nlowtot(t);
    Nreshigh(t) = freslow(t)*Nhightot(t);
    Nsenslow(t) = Nlowtot(t)-Nreslow(t);
    Nsenshigh(t) = Nhightot(t)-Nreshigh(t);
    
end

errorbarsnumsens = Nsenshigh-Nsenslow
errorbarsnumres = Nreshigh-Nreslow
errorbarsnumtot = Nhightot-Nlowtot



% Make into one column figure for single treatment
figure (3)
hold off
subplot(1,3,2)
hold off
errorbar(0:1:8, prolif_mean(1:9),error_bar_single_dose(2:10)/2, '-ro', 'LineWidth',2);
xlabel ('Time (weeks)', 'FontSize', 10)
ylabel ('Net Growth Rate per Day', 'FontSize', 18)
title('b', 'FontSize',18)
set(gca,'LineWidth',1.2,'FontSize',10)
%ylim([-0.5,0.5])
xlim([0 8])

subplot(1,3,1)
hold off
errorbar(0:1:8, fres2, errorbars_single_treat2_8./2, 'bo-', 'LineWidth', 1.2)
hold on
errorbar(0:1:8, fsens2, errorbars_single_treat2_8./2, 'yo-', 'LineWidth', 1.2)
xlabel ('Time (weeks)', 'FontSize', 10)
ylabel('Fraction of Cells','FontSize', 10)
title('a','FontSize', 18)
set(gca,'LineWidth',1.2,'FontSize',10)
xlim([0 8])
legend('Resistant', 'Sensitive')
legend boxoff

subplot(1,3,3)
hold off
errorbar(0:1:8, Nres_wks, errorbarsnumres, 'bo-', 'LineWidth',1.2)
hold on
errorbar(0:1:8, Nsens_wks, errorbarsnumsens, 'yo-', 'LineWidth',1.2)
%semilogy(0:1:8, Nhightot, 'ro-', 0:1:8, Nlowtot, 'ko-')
% semilogy(0:1:8, Ntot_wks,'g*-', 0:1:8, Nres_wks, 'b*-', 0:1:8, Nsens_wks, 'y*-', 'LineWidth', 2)
xlabel ('Time (weeks)', 'FontSize', 10)
ylabel ('Number of Cells', 'FontSize', 10)
set(gca,'LineWidth',1.2,'FontSize',10)
set(gca, 'Yscale', 'log')
xlim([0 8])
title('c','FontSize', 18)
legend('Resistant', 'Sensitive')
legend 'boxoff'
%%
figure(5)
hold off
subplot(1,3,1)
hold off
errorbar(0:1:8, prolif_mean(1:9),error_bar_single_dose(2:10)/2, '-ko', 'LineWidth',1.2);
xlabel ('Time (weeks)', 'FontSize', 10)
ylabel ('Per Capita Growth (cells per day)', 'FontSize', 18)
title('a', 'FontSize',18)
set(gca,'LineWidth',1.2,'FontSize',10)
%ylim([-0.5,0.5])
xlim([-0.1 8.1])
ylim([0.5 1.8])
legend boxoff

subplot(1,3,2)
hold off


hold on
errorbar(0:1:8, Nres_wks, errorbarsnumres, 'ro-', 'LineWidth',1.2)
errorbar(0:1:8, Nsens_wks, errorbarsnumsens, 'bo-', 'LineWidth',1.2)
%semilogy(0:1:8, Nhightot, 'ro-', 0:1:8, Nlowtot, 'ko-')
% semilogy(0:1:8, Ntot_wks,'g*-', 0:1:8, Nres_wks, 'b*-', 0:1:8, Nsens_wks, 'y*-', 'LineWidth', 2)
xlabel ('Time (weeks)', 'FontSize', 10)
ylabel ('Number of Cells', 'FontSize', 10)
set(gca,'LineWidth',1.2,'FontSize',10)
set(gca, 'Yscale', 'log')
xlim([-0.1 8.1])
ylim([10 1e11])
title('b','FontSize', 18)
legend('resistant', 'sensitive')
legend 'boxoff'
% Make insert
subplot(1,3,3)
hold off
errorbar(0:1:3, Nres_wks(1:4), errorbarsnumres(1:4), 'ro-', 'LineWidth',1.5)
hold on
errorbar(0:1:3, Nsens_wks(1:4), errorbarsnumsens(1:4), 'bo-', 'LineWidth',1.5)
xlabel ('Time (weeks)', 'FontSize', 10)
ylabel ('Number of Cells', 'FontSize', 10)
set(gca,'LineWidth',1.2,'FontSize',10)
% set(gca, 'Yscale', 'log')
xlim([0 3])
title('c','FontSize', 18)
%legend('Resistant', 'Sensitive')
legend 'boxoff'
