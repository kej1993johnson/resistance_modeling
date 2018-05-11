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
kobs = load('kobs_singletreat.mat'); % note will ignore first value since is before t = 0;
kobs = struct2cell(kobs);
kobs = cell2mat(kobs);
%% This part is specific to interval 1
Nin = [Nsens(1), Nres(1)];
Nin_tot = sum(Nin);
Nout = [Nsens(2),Nres(2)];
kobs = kobs(2:end);
% function [k pi] = Markov_transition(Nin, Nout, kobs)
% Given an initial resistant and sensitive cell number and a combined
% proliferation rate try to come up with solution to transition and
% individual ks

% I think your "energy" is given by the difference in Nout and the Nattempt

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
dt = 0.1;
t_steps = 7./dt
k_wk= (kobs(1)).^7  % This is the total proliferation in a week (so daily to the ^7)
k_mark = (log(k_wk))./7 % this should give you cells per day
t_int = 7;% interval is 7 days

% check that they are the same
Nout_tot_prolif = Nin_tot*exp(k_mark*t_int)
Nout_tot = sum(Nout) 


pisr1 = 0.8; % must be less than 1
pirs1 = 0.3; % must be less than 1
%%

% Sets up transition matrix
pi_init_guess = 0.5*ones(2,2);
pi_init_guess(1,2) = pisr1; % pi s->r
pi_init_guess(1,1) = 1-pi_init_guess(1,2); % pi s->s ****
pi_init_guess(2,1) = pirs1; % pi r->s
pi_init_guess(2,2) = 1-pi_init_guess(2,1); % pi r->r

% diagonals must add to kobs
% assume that sensitive cells are dying while resistant cells are able to
% proliferate
% we have the equation e^kdst + e^kdrt = e^kt

% Unconstrained k values ( arbitrary and not related to k
k_init_guess = zeros(2,2); 
kds1 = k_mark;
kds2 = k_mark;
k_init_guess(1,1) = kds1 ; % LHS diagonal kds
%k_init_guess(2,2) = (sum(Nin)./Nin(2))*(k_wk-k_init_guess(1,1)) + k_init_guess(1,1); % kdr scaled properly
k_init_guess(2,2) = kds2;
% Note that there are only three true unknowns here: pi s->r,pi r->s, and kds
% Therefore these are the paramters we need to optimize.
%%
N_guess = zeros(t_steps, 2);
N_guess(1,:) = Nin;



% Allow for state transitions
for t = 2:t_steps
    N_guess_curr = N_guess(t-1,:)';
    N_guess(t,:) = N_guess_curr + dt*(k_init_guess*pi_init_guess*(N_guess_curr));
    N_guess(t,:) = N_guess(t,:)';
end

% Do not allow for state transitions
k_nt = zeros(2,2);
N_guess_nt(1,:) = Nin;

for t = 2:t_steps
    kdsnt1= .001;
    k_nt(1,1) = kdsnt1;
    k_nt(2,2) = ((exp(k_mark*t*dt) - kdsnt1*t*dt)./(t*dt));
    N_guess_curr_nt = N_guess_nt(t-1,:)';
    N_guess_nt(t,:) = N_guess_curr_nt + dt*(k_nt*(N_guess_curr_nt));
    N_guess_nt(t,:) = N_guess_nt(t,:)';
end



Noutobs = sum(Nout)
Nout_trans = sum(N_guess_curr)
Nout_nt = sum(N_guess_curr_nt)
%%
figure(1)
plot(1:1:t_steps, log(N_guess(:,1)),'y')
hold on
plot(1:1:t_steps,log(N_guess(:,2)), 'b')
hold on
plot(1:1:t_steps,log(N_guess_nt(:,1)),'y.')
plot(1:1:t_steps, log(N_guess_nt(:,2)), 'b.')
plot(t_steps, log(Nout(1)), 'yo')
plot(t_steps, log(Nout(2)), 'bo')
plot(1, log(Nin(1)),'yo')
plot(1, log(Nin(2)),'bo')
