close all, clear all, clc
%%

fres2 = load('fres2.mat');
fres2 = struct2cell(fres2);
fres2 = cell2mat(fres2);
prolif = load('prolif_single_treatment.m'); % in Number of cells per day
% prolif = struct2cell(prolif);
% prolif = cell2mat(prolif);
prolif = vertcat(prolif, prolif(end,1));

prolif_reps = load('prolif_single_treatment_reps.m');

CIsingledose_8 =  prctile(prolif_reps', [2.5 97.5]);
error_bar_single_dose = CIsingledose_8(2,:) - CIsingledose_8(1,:);

k = prolif-1;
dt = .1;
t_wks = 9;
Ntot(1) = 1e6;
Nres(1) = Ntot(1)*fres2(1);
Nsens(1) = Ntot(1) - Nres(1);
num_steps = (t_wks*7)./dt;
%% FD so with k not constant. Need to resize k to be length of the simulation domain
k_large = imresize(k(2:10),[num_steps, 1],'nearest');
f_res_large = imresize(fres2,[num_steps, 1],'nearest');
%%
for t = 2:num_steps
    Ntot(t) = Ntot(t-1) + dt.*(k_large(t-1).*Ntot(t-1));
    Nres(t) = Ntot(t)*f_res_large(t);
    Nsens(t) = Ntot(t) - Nres(t);
    
end