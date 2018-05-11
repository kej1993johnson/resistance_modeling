% This script will be used to make a function that will find the optimal
% transition matrix and proliferation
close all, clear all, clc
% First load all 9 weeks of data, only will find transition matrix for one
% interval at a time though
Nres = load('Nres_wks_singletreat.mat');
Nres = struct2cell(Nres);
Nres = cell2mat(Nres);
Nsens = load('Nsens_wks_singletreat.mat');
Nsens = struct2cell(Nsens);
Nsens = cell2mat(Nsens);
kobs = load('kobs_singletreat.mat'); % note will ignore first value since is before t = 0
kobs = struct2cell(kobs); % also note, this is in terms of cells per day, need to subtract by 1 to get growth per day
kobs = cell2mat(kobs);
kobs = kobs-1; % daily proliferation rate since data obtained is in factor to multiply cells per day by!
fres = load('fres_wks_singletreat.mat');
fres = struct2cell(fres);
fres = cell2mat(fres);
fsens = 1-fres;
% This part is specific to interval 1
Nin = [Nsens(1), Nres(1)];
Nin_tot = sum(Nin);
Nout = [Nsens(2),Nres(2)];
Nout_tot = sum(Nout);
kobs = kobs(2:end);
% function [k pi] = Markov_transition(Nin, Nout, kobs)
% Given an initial resistant and sensitive cell number and a combined
% proliferation rate try to come up with solution to transition and
% individual ks

% I think your "energy" is given by the difference in Nout and the Nattempt
%% Try just using fres and fsens to find delta E
R = 8.32;
T = 278;
for t = 1:9
    deltaE(t) = log(fres(t)./fsens(t));
end

plot(0:1:8, deltaE, 'LineWidth', 3);
xlabel ('Time Post Treatment (Weeks)', 'FontSize', 10)
ylabel ('Es-Er (delta E)', 'FontSize', 10)
set(gca,'LineWidth',1.2,'FontSize',8)
title('Dynamic Change in Energy Difference Between States','FontSize', 10)

%%
% make a pi guess and a k guess based on inputs
% all rows add to 1
% Assume there is a strong high transition rate towards r and a low
% transition rate towards s

%First just set up an initial guess that correctly outputs the total cell number
% expect total cell number to be accurate because we know the initial total
% cell number and the the bulk observed k

% First step: write a function that finds Nguess
% kds must be limited to kobs^7 (for first time point this is 1.67)

%[Nguess_out] = find_Nguess(kds, pisr, pirs, Nin, kobs)
kobs1= kobs(1);
t_int = 7;% interval is 7 days
kds1 = 1.2*kobs1; % limit this (somewhat arbitrary?)
kdr1 = 0.8*kobs1;
pisr1 = 0.9; % must be less than 1
pirs1 = 0.2; % must be less than 1


% Check k mark is correct
%Nout_tot_prolif = Nin_tot*exp(kobs1*t_int);
Nout_tot = sum(Nout);
% Call function to find initial guess
Nout_guess1 = find_Nguess(kds1, kdr1, pisr1, pirs1, Nin, kobs1);
Nout_guess1tot = sum(Nout_guess1);

% Now it's just parameter optimization, potentially we can solve this
% analytically but we want to essentially find pi s->r,pi r->s, and kds so
% that we can minimize the difference between Nout_guess and Nout

% Find error in first guess
Err_init = sum((Nout_guess1-Nout).^2)./Nout_tot; 
Err_current = Err_init;
count_acc = 0;

%%
% start with this outside loop, this is relatively random to start with
% although ideally we get good at guessing this for each interval
N = 400;
T= 150;
%k = 1:1:N;
%T = T0* exp(-k./N*5);

%plot(k, T)

lambda = 0.01; % step size for searching pi
beta = 0.05*lambda; % step size for searching k
alpha = 5; % factor that kds and kdr can grow by
for k = 1:N

% The next chunk of code is the equivalent of setting Y = X and then
% changing Y and finding its new error
kds2 = kds1 + beta*(2*rand-1);
kdr2 = kdr1 + beta*(2*rand-1);% increments random number between -.1 and .1 ( note the .1 is arbitrary)
% Consider putting bounds on kds?
if kds2 > alpha*kobs1
    kds2 = alpha*kobs1;
end
if kds2 < -alpha*kobs1
    kds2 = -alpha*kobs1;
end

if kdr2 > alpha*kobs1
    kdr2 = alpha*kobs1;
end
if kdr2 < -alpha*kobs1
    kdr2 = -alpha*kobs1;
end

pisr2 = pisr1 + lambda*(2*rand-1);
if pisr2 >1
    pisr2 = 1;
end
if pisr2 <0
    pisr2 = 0;
end

pirs2 = pirs1 + lambda*(2*rand-1);
if pirs2 >1
    pirs2 = 1;
end
if pirs2 <0
    pirs2 = 0;
end
% want to store all pirs, kds, and psrs and error
kdse(k) = kds2;
kdre(k) = kdr2;
pisre(k) = pisr2;
pirse(k) =pirs2;

Nout_guess2 = find_Nguess(kds2,kdr2, pisr2, pirs2, Nin, kobs1);
Nout_guesse(k,1:2) = Nout_guess2(1:2);
%Err_new = (((Nout_guess2(1)-Nout(1)).^2) + ((Nout_guess2(2)-Nout(2)).^2))./Nout_tot; 
Err_new = sum((Nout_guess2-Nout).^2)./Nout_tot; 
Err_newe(k) = Err_new;
prob(k) = exp(((Err_current-Err_new))./T);
%prob_try(k) = exp(((Nout(1)-Nout_guess2(1)))./T);
%val(k) = 1+(2*rand-1)./T;
  if rand < exp(((Err_current-Err_new))./T) %if true, accept the change!
      % rand is any number between 0 and 1, and so this is tending to
      % always be
       kds1 = kds2;
       pirs1 = pirs2;
       pisr1 = pisr2;
       Nout_guess1 = Nout_guess2;
       Err_current = Err_new;
       count_acc= count_acc +1;
       lambda = 0.999*lambda;
       beta = 0.999*beta;
  else 
      lambda = 1.01*lambda;
      beta = 1.01*beta;

  end
    T = T*exp(-k./N*5);
    Te(k)=T;
    Score(k)=Err_current;
    kds1_acc(k) = kds1;
    kdr1_acc(k)= kdr1;
    pisr1_acc(k)=pisr1;
    pirs1_acc(k)=pirs1;
    Nout_guess_acc_all(k, :) = Nout_guess1;
end
%%

    figure(2)
    set(gca,'LineWidth',1.5,'FontSize',12);
    plot(1:1:N, kdse,'r*')
    hold on
    plot(1:1:N, kdre,'g*')
    plot(1:1:N, pisre,'b.')
    plot(1:1:N, pirse,'y.')
    plot(1:1:N,(Err_newe), 'k.')
    xlabel('Markov step')
 
    ylabel('Value')
    legend('kds', 'kdr','transition s->r', 'transition r->s', 'Error in prediction')
    title('Search of Parameter Space')
    
    % Next want to see kds, pisr, and pirs as a function of the error
    figure(3)
    set(gca,'LineWidth',1.5,'FontSize',12);
    plot(Err_newe, kdse, 'ro', 'LineWidth',2)
    hold on
    plot(Err_newe, kdre, 'go', 'LineWidth',2)
    plot(Err_newe, pisre, 'bo','LineWidth',2)
    plot(Err_newe, pirse, 'yo','LineWidth',2)
    xlabel('Error in prediction')
    ylabel('Parameter Value')
    legend ('kds','kdr', 'transition s->r', 'transition r->s')
    title('Parameter value as a function of error')
    %xlim([0, 0.03])
    %ylim([0, 1.5])
    
    figure(4)
    set(gca,'LineWidth',1.5,'FontSize',12);
    plot(kdse, Err_newe, 'r*')
    hold on
    plot(kdre, Err_newe, 'g*')
    plot(pisre, Err_newe, 'b.')
    plot(pirse, Err_newe, 'y.')
    xlabel('Parameter value')
    ylabel('Error in prediction')
    legend ('kds', 'kdr','transition s->r', 'transition r->s')
    title('Error as a function of parameter value')
    
    
    Nout_all = Nout.*ones(N,2);
    figure(6)
    set(gca,'LineWidth',1.5,'FontSize',12);
    plot(1:1:N, Nout_guesse(:,1), 'y.')
    hold on
    plot (1:1:N, Nout_guesse(:,2), 'b.')
    plot(1:1:N, Nout_all(:,1), '-y')
    plot(1:1:N, Nout_all(:,2),'-b')
    
    figure(7)
    set(gca,'LineWidth',1.5,'FontSize',12);
    plot(1:1:N,prob, 'b.')
    hold on 
    title('Distribution of Probabilities')
   % plot(1:1:N, val, 'r.')
    
    figure(8)
    set(gca,'LineWidth',1.5,'FontSize',12);
    plot(1:1:N, kds1_acc, 'r.')
    hold on
    plot(1:1:N, kdr1_acc, 'g.')
    plot(1:1:N, pisr1_acc,'b.')
    plot(1:1:N, pirs1_acc,'y.')
    legend('gs', 'gr', 'transition s->r', 'transition r->s')
    
   figure(6)
    set(gca,'LineWidth',1.5,'FontSize',12);
    plot(1:1:N, Nout_guess_acc_all(:,1), 'y.')
    hold on
    plot (1:1:N, Nout_guess_acc_all(:,2), 'b.')
    plot(1:1:N, Nout_all(:,1), '-y')
    plot(1:1:N, Nout_all(:,2),'-b')
    
 
    
%     figure(6)
%     hold off
%     surf(pisre, pirse, Z);
    






