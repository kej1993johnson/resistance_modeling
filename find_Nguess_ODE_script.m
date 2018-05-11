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
% function [k pi] = Markov_transition(Nin, Nout, kobs)
% Given an initial resistant and sensitive cell number and a combined
% proliferation rate try to come up with solution to transition and
% individual ks

k_wk= (kobs(1)).^7; % This is the total proliferation in a week (so daily to the ^7)
k_mark = (log(k_wk))./7 ;% this should give you cells per day
t_int = 7;% interval is 7 days

% initial guesses to be put into function
gs = 1.5*k_mark;
gr = 0.4*k_mark;
krs = 0.9;
ksr = 1.0;

dt = 0.001;
t_steps = 7./dt;


% Set up growth matrix
G = zeros(2,2);
G(1,1) = gs;
G(2,2) = gr;

% Set up transition rate matrix
T = zeros(2,2);
T(1,1) = -ksr;
T(2,1) = ksr; % krr
T(2,2) = -krs;
T(1,2) = krs; % kss


N_guess = zeros(t_steps, 2);
N_guess(1,:) = Nin;



% Allow for state transitions
for t = 2:t_steps
    N_guess_curr = N_guess(t-1,:)';
    N_guess(t,:) = N_guess_curr + dt*((G + T)*(N_guess_curr));
end

Nout_guess_fin = N_guess(t_steps, :); % outputs final guess 

plot(1:1:t_steps, N_guess(:,1), 'y')
hold on
plot(1:1:t_steps, N_guess(:,2), 'b')