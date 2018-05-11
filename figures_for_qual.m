clear all, close all, clc
% This script makes figures for qual:

% First make scatterplot of large variance in x and y

ax = 5;
bx1 = 35;
ay = 5;
by1 = 35;
bx2 = 60;
by2 = 60;
y1 = ay.*randn(500,1) + by1;
x1 = ax.*randn(500,1) + bx1;
y2 = ay.*randn(500,1) + by2;
x2 = ax.*randn(500,1) + bx2;


ax = 10;
bx1 = 35;
ay = 10;
by1 = 35;
bx2 = 60;
by2 = 60;
y5 = ay.*randn(500,1) + by1;
x5 = ax.*randn(500,1) + bx1;
y6 = ay.*randn(500,1) + by2;
x6 = ax.*randn(500,1) + bx2;

figure(1)
subplot(1,2,1)
plot(x1,y1, 'b.')
hold on
plot(x2,y2, 'b.')
set(gca,'LineWidth',1.5,'FontSize',12);
xlabel('PC1')
ylabel('PC2')
xlim([0 100])
ylim([0 100])
title('Two Low-Variability Subpopulations')
subplot(1,2,2)
plot(x5,y5, 'b.')
hold on
plot(x6,y6, 'b.')
set(gca,'LineWidth',1.5,'FontSize',12);
xlabel('PC1')
ylabel('PC2')
xlim([0 100])
ylim([0 100])
title('Two High-Variability Subpopulations')


%%
ax3 = 15;
bx3 = 50;
by3 = 50;
ay3 = 15;
x3 = ax3.*randn(500,1) + bx3;
y3 = ay3.*randn(500,1) + by3;
ax4 = 5;
ay4 = 5;
x4 = ax4.*randn(500,1) + bx3;
y4 = ay4.*randn(500,1) + by3;

figure(2)
subplot(1,2,1)
plot(x3,y3, 'r.')
xlim([0 100])
ylim([0 100])
xlabel('PC1')
ylabel('PC2')
title('High Variability in Gene Expression')
subplot(1,2,2)
plot(x4,y4, 'r.')
xlim([0 100])
ylim([0 100])
xlabel('PC1')
ylabel('PC2')
title('Well-Defined Gene Expression')
%%
x = -15:.1:15;
a1 = -50;
a2 = -90;
F1 = x.^4 + 3.*(x.^3) +a1.*(x.^2) ;
F2 = x.^4 -1.*(x.^3) +a2.*(x.^2) ;

figure(3)
hold off
set(gca,'LineWidth',1.5,'FontSize',12);
set(gca,'ytick',[])
set(gca,'xtick',[])
plot(x, F1, 'b-','LineWidth', 3)
hold on
plot(x, F2, 'y-', 'LineWidth', 3)
xlabel('Drug Sensitivity')
ylabel('Free Energy')
ylim([-3000, 1000])
xlim([-10, 12])
