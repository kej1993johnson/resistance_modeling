subplot(1,2,1) % graphs fit to all weeks in red
hold off
hold on
plot (X, v_model1plot, '-r', 'LineWidth', 3);
plot(X, v_model1plotonedose, '-g', 'LineWidth',3);
plot(X, v_model1plotunt, '-k', 'LineWidth',3);
xlim([0,250]);
xlabel('Dose (uM)','FontSize',18)
ylabel('Viability','FontSize',18)
title ('Single Static Population Model Effect of Pulse Treatment','FontSize',14)
legend('Twice Treated LD50 = 57.47 +/- 3.6 uM', 'Single Treatment LD50 = 51.30 +/- 1.75 uM', 'Untreated LD50 = 38.18 +/- 1.2 uM')
set(gca,'LineWidth',1.5,'FontSize',10)

time = 0:1:15;
subplot(1,2,2)
hold off
% Plot of LD50s over time
%plot(time, betaLD50(16:30), 'o', 'LineWidth', 4)
errorbar(time, fres2twodoseall, errorbar2twodose./2, '-ro', 'LineWidth', 2)
hold on
errorbar(time,fres2onedoseall, errorbar2onedose./2, '-go', 'LineWidth',2)
xlabel ('Time Post Treatment (Weeks)', 'FontSize', 14)
ylim([0 1])
ylabel ('Fraction of Cells Resistant', 'FontSize', 14)
set(gca,'LineWidth',1.5,'FontSize',10)
%errorbar(time, fres2untall, errorbar2untreated./2,'ko', 'LineWidth',2)
title('Two Population Model Effect of Pulse Treatments ','FontSize', 12)
legend('Twice Treated ','Single Treatment', 'Untreated Resistant Fraction')
