% Script finds and plots the average intensity in the cytoplasm for each
% ploidy cell type over time from the incucyte

close all, clear all ,clc
ploidy_data = load('ploidy_data.m');
% data is arranged so that top 8 rows on left ( 5 columns) are total
% integrated intensity per well with low threshold, right of that is the
% total green object area, then below is the same for the high threshold

tot_cell_int = ploidy_data(1:8, 1:5);
total_cell_area = ploidy_data(1:8, 6:10);
nuc_cell_int = ploidy_data(9:16, 1:5);
nuc_cell_area = ploidy_data(9:16, 6:10);

cyt_cell_area = total_cell_area- nuc_cell_area;
cyt_int = tot_cell_int-nuc_cell_int;

avg_int_cyt = cyt_int./cyt_cell_area;

t = 0:0.25:1.75;

whole_pop = avg_int_cyt(:,1);
cyt_int_8N = avg_int_cyt(:,2);
cyt_int_4N = avg_int_cyt(:,3);
cyt_int_2N = avg_int_cyt(:,4);
untreated = avg_int_cyt(:,5);

figure(1)
xlim([0,1.75]);
plot(t, whole_pop, t, cyt_int_8N, t, cyt_int_4N, t, cyt_int_2N, 'LineWidth', 2)
xlabel('Time after dye injection (hours)','FontSize',18)
ylabel('Average Dye Intensity in Cytoplasm','FontSize',18)
title ('Ploidy Sorted Dye Pumping Assay','FontSize',18)
legend('Whole Population','8N','4N', '2N')
set(gca,'LineWidth',1.5,'FontSize',12)

final_time = avg_int_cyt(5,1:4);
final_area = total_cell_area(5,1:4);
final_int = cyt_int(5,1:4);

figure(2)
hold off
bar(final_time)
xlabel('Ploidy Sorted Populations','FontSize',18)
ylabel('Average Cytoplasmic Dye Intensity at Final Time','FontSize',18)
title('Dye Intensity at Final Time','FontSize',18)
set(gca,'LineWidth',1.5,'FontSize',12)

figure(3)
hold off
bar(final_area)
xlabel('ploidy sorted')
ylabel('total cell area at final time')
set(gca,'LineWidth',1.5,'FontSize',12)

figure(3)
hold off
bar(final_int)
xlabel('ploidy sorted')
ylabel('total cell int at final time')
set(gca,'LineWidth',1.5,'FontSize',12)

