% This script will be used to make a function that will find the growth and
% transition rates using Monte Carlo methods to search parameter space and
% ODE based model
close all, clear all, clc
% First load all 9 weeks of data, only will find transition matrix for one
% interval at a time though
Nres = load('Nres_wks_singletreat.mat');
Nres = struct2cell(Nres);
Nres = cell2mat(Nres);
Nsens = load('Nsens_wks_singletreat.mat');
Nsens = struct2cell(Nsens);
Nsens = cell2mat(Nsens);
kobs = load('kobs_singletreat.mat'); % note will ignore first value since is before t = 0;
kobs = struct2cell(kobs);
kobs = cell2mat(kobs);
% This part is specific to interval 1
Nin = [Nsens(1), Nres(1)];
Nin_tot = sum(Nin);
ratio_in = Nin(1)./Nin(2);
Nout = [Nsens(2),Nres(2)];
Nout_tot = sum(Nout);
ratio_out = Nout(1)./Nout(2);
kobs = kobs(2:end);
k_exp = kobs-1;
% function [k pi] = Markov_transition(Nin, Nout, kobs)
% Given an initial resistant and sensitive cell number and a combined
% proliferation rate try to come up with solution to transition and
% individual ks

k_wk= (kobs(1)).^7; % This is the total proliferation in a week (so daily to the ^7)
k_mark = (log(k_wk))./7 ;% this should give you cells per day
t_int = 7;% interval is 7 days
k_mark = kobs(1)-1;
% initial guesses to be put into function
gs1 = 1.1*k_mark;
gr1 = 0.9*k_mark;
krs1 = 0.5;
ksr1 = 0.55;

Nguess_ODE = find_Nguess_ODE(gs1,gr1,ksr1,krs1, Nin, k_mark);
Nguess_tot = sum(Nguess_ODE);
ratio_guess = Nguess_ODE(1)./Nguess_ODE(2);
Nout_guess1 = Nguess_ODE;

err_ct_init =abs(( Nout_tot - Nguess_tot)./Nout_tot);
err_ratio_init = abs(ratio_out-ratio_guess);
err_ct_curr = err_ct_init;
err_ratio_curr = err_ratio_init;

err_all_init = (abs((Nout(1)-Nguess_ODE(1))./Nout(1)) + abs((Nout(2)-Nguess_ODE(2))./Nout(2)));
err_all_curr = err_all_init;
count_acc = 0;
% Find that the initial guess has almost equal parts error in it's ratio
% and it's count. Now we want to search parameter space to change the
% unknowns gs, gr, krs, and ksr to decrease both of these errors
N = 2000;
n_wks = 8;
store_params = zeros(N, 7, n_wks);

store_acc_val = zeros(n_wks, 7);
%%

Gs=-log(Nsens./(Nres+Nsens));
Gr = -log(Nres./(Nres+Nsens));
deltaG = -log(Nres./Nsens);
hold off
plot(0:1:8, Gs, 'y')
hold on
plot(0:1:8, Gr, 'b')
plot(0:1:8, deltaG, 'k')

%%  Perform Monte Carlo method to search parameter space

% start with this outside loop, this is relatively random to start with
% although ideally we get good at guessing this for each interval

T0 = 0.2;
k=1:1:N;
T = T0*exp(-k./N*5);% sets up gradual "cooling"

figure(8)
hold off
set(gca,'LineWidth',1.5,'FontSize',12);
plot(k, T, 'LineWidth', 3)
xlabel('Markov step')
ylabel('Temperature')
title('Temperature Annealing')


lambda = 0.03; % step size for searching pi
beta = lambda; % step size for searching k
alpha = 3; % factor that kds and kdr can grow by *** This is important, we need to ask about this.
           % What is physiological variance in growth rates of different
           % cell types??

%%
for k = 1:N

% The next chunk of code is the equivalent of setting Y = X and then
% changing Y and finding its new error
gs2 = gs1 + beta*(2*rand-1);
gr2 = gr1 + beta*(2*rand-1);% increments random number between -.1 and .1 ( note the .1 is arbitrary)
% Consider putting bounds on kds?


if gs2 > alpha*k_mark
    gs2 = alpha*k_mark;
end

if gs2 < -alpha*k_mark
     gs2 = -alpha*k_mark;
end

if gr2 > alpha*k_mark
    gr2 = alpha*k_mark;
end

if gr2 < -alpha*k_mark
     gr2 = -alpha*k_mark;
end

if gr2 <0
    if gr2<gs2
        gr2= gs2;
    end
end


    

ksr2 = ksr1 + lambda*(2*rand-1);
if ksr2 >1
    ksr2 = 1;
end
if ksr2 <0
    ksr2 = 0;
end

krs2 = krs1 + lambda*(2*rand-1);
if krs2 >1
    krs2 = 1;
end
if krs2 <0
    krs2 = 0;
end
% want to store all pirs, kds, and psrs and error
gse(k) = gs2;
gre(k) = gr2;
ksre(k) = ksr2;
krse(k) =krs2;

Nout_guess2 = find_Nguess_ODE(gs2,gr2,ksr2,krs2, Nin, k_mark);
Nguess_tot2= sum(Nout_guess2);
Nout_guesse(k,1:2) = Nout_guess2(1:2);
ratio_guess2 = Nout_guess2(1)./Nout_guess2(2);
err_ct_new =abs(( Nout_tot - Nguess_tot2)./Nout_tot);
err_ratio_new = abs(ratio_out-ratio_guess2);
err_all_new = (abs((Nout(1)-Nout_guess2(1))./Nout(1)) + abs((Nout(2)-Nout_guess2(2))./Nout(2)));
err_all_newe(k) = err_all_new;

prob(k) = exp(((err_all_curr-err_all_new))./T(k));
%prob_try(k) = exp(((Nout(1)-Nout_guess2(1)))./T);
%val(k) = 1+(2*rand-1)./T;

  if rand < prob(k) %if true, accept the change!
      % if new error is less than current error, then rand between 0 and 1
      % always less than prob since prob will have a positive number in the
      % exponent so it will be greater than 1
      % if prob has a negative number in the exponent (so current error is less than new, we want to sometimes
      % accept but not usually
       gs1 = gs2;
       gr1 = gr2;
       krs1 = krs2;
       ksr1 = ksr2;
       Nout_guess1= Nout_guess2;
       err_ratio_curr = err_ratio_new;
       err_ct_curr = err_ct_new;
       err_all_curr = err_all_new;
       count_acc= count_acc +1;
       lambda = 0.999*lambda;
       beta = 0.999*beta;
  else 
      lambda = 1.001*lambda;
      beta = 1.001*beta;

  end

    Score_ratio(k)=err_ratio_curr;
    Score_tot(k)= err_ct_curr;
    Score(k) = err_all_curr;
    gr_acc(k) = gr1;
    gs_acc(k) = gs1;
    krs_acc(k)= krs1;
    ksr_acc(k) =ksr1; 
    Nout_guess_acc(k,:) = Nout_guess1;
    
end
%%
store_params1_uns = (horzcat(ksre', krse', gse', gre', Nout_guess_acc(:,1), Nout_guess_acc(:,2), err_all_newe'));
store_params1_sort = sortrows(store_params1_uns, 7);
store_params(:,:,1) = store_params1_sort;
store_acc_val(1,:) = horzcat(ksr1, krs2, gs1, gr1, Nout_guess1(1), Nout_guess1(2), err_all_curr);


%%
figure(1)
hold off
set(gca,'LineWidth',1.5,'FontSize',12);
plot(1:1:N,prob, 'b.')
hold on 
title('Distribution of Probabilities')



    figure(2)
    hold off
    set(gca,'LineWidth',1.5,'FontSize',12);
    plot(1:1:N, gr_acc,'b*')
    hold on
    plot(1:1:N, gs_acc, 'y*')
    plot(1:1:N, ksr_acc,'b.')
    plot(1:1:N, krs_acc,'y.')
    hold on
    plot(1:1:N,5*Score, 'r.')
    xlabel('Markov step')
    ylabel('Value')
    legend('gr', 'gs','transition s->r', 'transition r->s', 'Total Error')
    title('Search of Parameter Space')
    
    Nout_true = ones(N,2).*Nout;
    
    figure(3)
    hold off
    set(gca,'LineWidth',1.5,'FontSize',12);
    plot(1:1:N, Nout_true(:,1), 'y')
    hold on
    plot(1:1:N, Nout_true(:,2),'b')
    plot(1:1:N, Nout_guess_acc(:,1), 'y*')
    plot(1:1:N, Nout_guess_acc(:,2), 'b*')
    ylim([8e5, 9e5])
    ylabel('Cell Number')
    xlabel('Monte Carlo Step')
%%    
    % Next want to see kds, pisr, and pirs as a function of the error
    figure(4)
    set(gca,'LineWidth',1.5,'FontSize',12);
    plot(err_all_newe, krse, 'r.', 'LineWidth',2)
    hold on
    plot(err_all_newe, ksre, 'g.', 'LineWidth',2)
    plot(err_all_newe, gre, 'b.','LineWidth',2)
    plot(err_all_newe, gse, 'y.','LineWidth',2)
    xlabel('Error in prediction')
    ylabel('Parameter Value')
    legend ('krs','ksr', 'gre', 'gse')
    title('Parameter value as a function of error')
    xlim([0, 0.2])
    ylim([0, 1.5])
    
    figure(5)
    set(gca,'LineWidth',1.5,'FontSize',12);
    plot(krse, err_all_newe, 'r.', 'LineWidth',2)
    hold on
    plot(ksre, err_all_newe, 'g.', 'LineWidth',2)
    plot(gre, err_all_newe, 'b.','LineWidth',2)
    plot(gse, err_all_newe, 'y.','LineWidth',2)
    ylabel('Error in prediction')
    xlabel('Parameter Value')
    legend ('krs','ksr', 'gre', 'gse')
    title('Error as a function of paramter values')
    xlim([-1, 2])
    ylim([0, 0.05])
  
  %%  
    figure(4)
    set(gca,'LineWidth',1.5,'FontSize',12);
    plot(gse, err_all_newe, 'r*')
    hold on
    plot(gre, err_all_newe, 'g*')
    plot(ksre, err_all_newe, 'b.')
    plot(krse, err_all_newe, 'y.')
    xlabel('Parameter value')
    ylabel('Error in prediction')
    legend ('kds', 'kdr','transition s->r', 'transition r->s')
    title('Error as a function of parameter value')
    ylim([0 0.2])
    
    
    Nout_all = Nout.*ones(N,2);
    figure(6)
    set(gca,'LineWidth',1.5,'FontSize',12);
    plot(1:1:N, Nout_guesse(:,1), 'y.')
    hold on
    plot (1:1:N, Nout_guesse(:,2), 'b.')
    plot(1:1:N, Nout_all(:,1), '-b')
    plot(1:1:N, Nout_all(:,2),'-y')
    
    figure(7)
    set(gca,'LineWidth',1.5,'FontSize',12);
    plot(1:1:N,prob, 'b.')
    hold on 
    title('Distribution of Probabilities')
   % plot(1:1:N, val, 'r.')

%% Now want to write a code that does this for all 8 time intervals.
% Loop needs to store the top values of ks and gs, their model outcome, and
% the error associated with them

%reset Nout to a two column vector of Nres and Nsens
for i = 2:length(Nres)
    Ntot_all(i-1,1) = Nres(i) + Nsens(i);
end
Nout_all = horzcat(Nres(2:9),Nsens(2:9), Ntot_all);

%reset kmark
k_exp = kobs-1;
count_acc_all = zeros(8,1);
prob = zeros([N, 8]);
%%
% First store all ksr, krs, gs, gr, Nres, Nsens, and error in a matrix at end of loop
n_wks = 8;

store_params = zeros(N, 7, n_wks);

for n = 1:n_wks % loop for new time interval
    
    % Write out whole code from above but at end store parameters, and set 
    % the initial conditions of the next time = to the starting conditions
    % of the first.
    
    for k = 1:N % loop for each new search of parameter space

        % The next chunk of code is the equivalent of setting Y = X and then
        % changing Y and finding its new error
        gs2 = gs1 + beta*(2*rand-1);
        gr2 = gr1 + beta*(2*rand-1);% increments random number between -.1 and .1 ( note the .1 is arbitrary)
        % Consider putting bounds on kds?
        if gs2 > alpha*k_exp(n)
        gs2 = alpha*k_exp(n);
        end
        if gs2 < -alpha*k_exp(n)
        gs2 = -alpha*k_exp(n);
        end

        if gr2 > alpha*k_exp(n)
        gr2 = alpha*k_exp;
        end
        if gr2 < -alpha*k_exp(n)
        gr2 = -alpha*k_exp(n);
        end

    ksr2 = ksr1 + lambda*(2*rand-1);
       if ksr2 >1
        ksr2 = 1;
       end
       if ksr2 <0
        ksr2 = 0;
       end

    krs2 = krs1 + lambda*(2*rand-1);
    if krs2 >1
        krs2 = 1;
    end
    if krs2 <0
        krs2 = 0;
    end
    % want to store all pirs, kds, and psrs and error
    gse(k,n) = gs2;
    gre(k,n) = gr2;
    ksre(k,n) = ksr2;
    krse(k,n) =krs2;

    Nout_guess2 = find_Nguess_ODE(gs2,gr2,ksr2,krs2, Nin, k_exp(n));
    Nout_guess2sens(k,n) = Nout_guess2(1);
    Nout_guess2res(k,n) = Nout_guess2(2);
    Nguess_tot2= sum(Nout_guess2);
    Nout_guesse(k,1:2) = Nout_guess2(1:2);
    ratio_guess2 = Nout_guess2(1)./Nout_guess2(2);

    err_all_new = (abs((Nout_all(n,1)-Nout_guess2(1))./Nout_all(n,3)) + abs((Nout_all(n,2)-Nout_guess2(2))./Nout_all(n,3)));
    err_all_newe(k,n) = err_all_new;

prob(k,n) = exp(((err_all_curr-err_all_new))./T(k));
%prob_try(k) = exp(((Nout(1)-Nout_guess2(1)))./T);
%val(k) = 1+(2*rand-1)./T;

  if rand < prob(k,n) %if true, accept the change!
      % if new error is less than current error, then rand between 0 and 1
      % always less than prob since prob will have a positive number in the
      % exponent so it will be greater than 1
      % if prob has a negative number in the exponent (so current error is less than new, we want to sometimes
      % accept but not usually
       gs1 = gs2;
       gr1 = gr2;
       krs1 = krs2;
       ksr1 = ksr2;
       Nout_guess1= Nout_guess2;
       err_all_curr = err_all_new;
       count_acc_all(n,1)= count_acc_all(n,1) +1;
       lambda = 0.999*lambda;
       beta = 0.999*beta;
  else 
      lambda = 1.001*lambda;
      beta = 1.001*beta;

  end

    Score_ratio(k,n)=err_ratio_curr;
    Score_tot(k,n)= err_ct_curr;
    Score(k,n) = err_all_curr;
    gr_acc(k,n) = gr1;
    gs_acc(k,n) = gs1;
    krs_acc(k,n)= krs1;
    ksr_acc(k,n) =ksr1; 
    Nout_guess_acc(k,:,n) = Nout_guess1;
    
    % Want all "2"s to be stored, at end of k loop sort by lowest error and
    % add to score_params
    
    
    end
%store all ksr, krs, gs, gr, Nres, Nsens, and error in a matrix at end of loop
% store_params(:,:,n) =  horzcat(ksre(:,n), krse(:,n), gse(:,n), gre(:,n), Nout_guess2sens(:,n), Nout_guess2res(:,n), err_all_newe(:,n));
 
% reset lambda, alpha and beta each iteration
lambda = 0.03; % step size for searching pi
beta = lambda; % step size for searching k
alpha = 5; 
% Reset the initial guesses so that they reflect the observed ks
% gs1 = 1*k_exp(n);
% gr1 = 0.9*k_exp(n);
% krs1 = 0.5;
% ksr1 = 0.75;
% 
% Nguess_ODE = find_Nguess_ODE(gs1,gr1,ksr1,krs1, Nin, k_exp(n));
% Nguess_tot = sum(Nguess_ODE);
% ratio_guess = Nguess_ODE(1)./Nguess_ODE(2);
% Nout_guess1 = Nguess_ODE;
% 
% err_all_init = (abs((Nout_all(n,1)-Nout_guess1(1))./Nout_all(n,3)) + abs((Nout_all(n,2)-Nout_guess1(2))./Nout_all(n,3)));
% err_all_curr = err_all_init;


    
end
%%
for k = 1:n
    store_params_sorted(:,:,k) = sortrows(store_params(:,:,n),7);
end

%plot(1:1:8, store_params_sorted(1, 1,:))
params_best = store_params_sorted(1,:,:);
for i = 1
    params_best_corr(:,:,i) = params_best(1,:,:);
end
