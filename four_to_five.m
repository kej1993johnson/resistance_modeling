% This script finds the parameters of best fit from interval of 4 to 5
% weeks.

Nin = [Nsens(5), Nres(5)]; % 4 wks
Nin_tot = sum(Nin);
ratio_in = Nin(1)./Nin(2); 
Nout = [Nsens(6),Nres(6)]; % 5 wks
Nout_tot = sum(Nout);
ratio_out = Nout(1)./Nout(2);

% function [k pi] = Markov_transition(Nin, Nout, kobs)
% Given an initial resistant and sensitive cell number and a combined
% proliferation rate try to come up with solution to transition and
% individual ks


k_mark = k_exp(5);
% initial guesses to be put into function
gs1 = 1.7*k_mark;
gr1 = 0.6*k_mark;
krs1 = 0.5;
ksr1 = 0.7;

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
%%  Perform Monte Carlo method to search parameter space

% start with this outside loop, this is relatively random to start with
% although ideally we get good at guessing this for each interval
N = 2000;
T0 = 0.5;
k=1:1:N;
T = T0*exp(-k./N*5);% sets up gradual "cooling"
hold off
plot(k, T)
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
if abs(gs2) > abs(alpha*k_mark)
    if gs2<0
    gs2 = alpha*k_mark;
    end
    if gs2>0
    gs2 = -alpha*k_mark;
    end
end

 if gr2<0
    if gr2 <  gs2
    gr2 = gs2;
    end
 end  
    if gr2>0
    if gr2>-0.5.*alpha*k_mark
    gr2 = -0.5*alpha*k_mark;
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
store_params5_uns = (horzcat(ksre', krse', gse', gre', Nout_guess_acc(:,1), Nout_guess_acc(:,2), err_all_newe'));
store_params5_sort = sortrows(store_params5_uns, 7);
store_params(:,:,5) = store_params5_sort;
store_acc_val(5,:) = horzcat(ksr1, krs2, gs1, gr1, Nout_guess1(1), Nout_guess1(2), err_all_curr);
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
    %ylim([8e5, 9e5])
%%    
    % Next want to see kds, pisr, and pirs as a function of the error
    figure(4)
    hold off
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
    ylim([-.5, 1.5])
    
    figure(5)
    hold off
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