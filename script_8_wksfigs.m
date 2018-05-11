
subplot(1,3,1) % graphs fit to all weeks in red
hold off
plot (dose, var, 'g.', 'LineWidth', 1.5)
hold on
plot(dose,var_unt, 'k.', 'LineWidth', 1.5)
plot (X, v_model1plot, '-g', 'LineWidth', 2);
plot(X, v_model1plotunt, '-k', 'LineWidth',2);
xlim([0,250]);
xlabel('Dose (uM)','FontSize',14)
ylabel('Viability','FontSize',14)
title ('Single Static Population Model','FontSize',14)
legend('Treated Viability','Untreated Viability','Treated LD50 = 50.4 uM', 'Untreated LD50 = 37.0 uM')
set(gca,'LineWidth',1.5,'FontSize',10)

subplot(1,3,2)
hold off
% Plot of LD50s over time
%plot(time, betaLD50(16:30), 'o', 'LineWidth', 4)
errorbar(time, LD50s, errorbarsLD50./2, '-go', 'LineWidth', 1.5)
hold on
%errorbar(0, betaLD50naive(2), errorbarlength1bootnaive(2)./2,'ko', 'LineWidth',2)
%plot(0, betaLD50naive(2), 'go', 'LineWidth',4)
xlabel ('Time Post Treatment(Weeks)', 'FontSize', 14)
ylabel ('Single Population LD50 (uM)', 'FontSize', 14)
set(gca,'LineWidth',1.5,'FontSize',10)
errorbar(time, LD50sunt, errorbarsLD50unt./2,'ko', 'LineWidth',1.5)
ylim ([20,90])
title('Single Dynamic Population Model','FontSize', 14)
legend('Treated', 'Untreated')

subplot(1,3,3)
hold off
errorbar(time, fres2, errorbars./2, 'bo-', 'LineWidth', 1.5)
hold on
errorbar(time, fsens2, errorbars./2, 'yo-', 'LineWidth', 1.5)
%errorbar(time, fres2unt, errorbars_unt./2, 'ko', 'LineWidth', 1.5)
%errorbar(0, fres0, errorbarlength2bootnaive./2, 'bo', 'LineWidth', 1.5)
%errorbar(time, fres2unt, errorbars_unt./2, 'ko-', 'LineWidth', 2)
%errorbar(0, fsens0, errorbarlength2bootnaive./2, 'ro', 'LineWidth', 1.5)
ylim([ 0 1])
xlabel('Time Post Treatment(Weeks)', 'FontSize', 14)
ylabel('Fraction of Cells', 'FontSize',14)
title ('Two Population Model', 'FontSize', 14)
set(gca,'LineWidth',1.5,'FontSize',10)
legend (' Resistant LD50 = 79.7 uM', 'Sensitive LD50 = 22.4 uM')
hold off

% subplot(2,2,4)
% hold off
% hold on
% errorbar(time, fsens2, errorbars./2, 'yo-', 'LineWidth', 1.5)
% %errorbar(time, fres2unt, errorbars_unt./2, 'ko', 'LineWidth', 1.5)
% %errorbar(0, fres0, errorbarlength2bootnaive./2, 'bo', 'LineWidth', 1.5)
% errorbar(time, fsens2unt, errorbars_unt./2, 'ko', 'LineWidth', 2)
% %errorbar(0, fsens0, errorbarlength2bootnaive./2, 'ro', 'LineWidth', 1.5)
% ylim([ 0 1])
% xlabel('Time Post Treatment (Weeks)', 'FontSize', 14)
% ylabel('Fraction of Cells', 'FontSize',14)
% title ('Two Population Model: Fraction of Sensitive Cells', 'FontSize', 14)
% set(gca,'LineWidth',1.5,'FontSize',10)
% legend ('Fraction of Treated Cells Sensitive LD50 = 22.4 uM', 'Fraction of Untreated Cells Sensitive')
% hold off
%%
figure(7)
subplot(1,2,1)% looks at example weeks 3, 6, and 10 WPT against all experimental dat
hold off
plot (dose(sum(nsize(1))+1:sum(nsize(1:2))), var(sum(nsize(1))+1:sum(nsize(1:2))), 'ro')
hold on
plot (dose(sum(nsize(1:7))+1:sum(nsize(1:8))), var(sum(nsize(1:7))+1:sum(nsize(1:8))), 'co')
%plot (dose(sum(nsize(1:9))+1:sum(nsize(1:10))), var(sum(nsize(1:9))+1:sum(nsize(1:10))), 'bo')
plot (X, v_model2_2, '-r', 'LineWidth', 3);
plot (X, v_model2_8, '-c', 'LineWidth', 3);
%plot (X, v_model2_10, '-b', 'LineWidth', 3);
xlim ([0 250])
set(gca,'LineWidth',1.5,'FontSize',10)
xlabel('Dose (uM)','FontSize',14)
ylabel('Viability','FontSize',14)
legend('Measured Viability at 2 weeks', 'Measured Viability at 8 weeks', 'Model Fit at Two Weeks', 'Model Fit at Eight Weeks')
title ('Two Population Model Example Fit at 2 and 8 Weeks Post Treatment','FontSize',12)
%
y = [fres2(3), fsens2(3); fres2(9), fsens2(9)];
z = [2;8];
err =[ errorbars(3), errorbars(3); errorbars(9), errorbars(9)];
c = categorical({ '2 weeks', '8 weeks'});
x(:,1) = fres2;
x(:,2) = fsens2;


subplot(1,2,2)
hold off
bar(z,y, 0.5, 'stacked')
hold on
set(gca,'LineWidth',1.5,'FontSize',10)
ylabel('Fraction of Population','FontSize',14)
xlabel(' Weeks Post Treatment', 'FontSize',14)
ylim ([0 1])
legend('Resistant Fraction', 'Sensitive Fraction')
title('Subpopulation Estimates at 2 and 8 Weeks','FontSize',12)
%%
b(:,1) = fres2;
b(:,2) = fsens2;
figure(9)
bar3( time, b, 'stacked')
set(gca,'LineWidth',1.5,'FontSize',10)
ylabel('Weeks Post Treatment', 'FontSize', 12)
xlabel('Subpopulation')
zlabel('Fraction of Cells')
legend('Fraction Sensitive', 'Fraction Resistant')

