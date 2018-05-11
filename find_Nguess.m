function [Nout_guess_fin] = find_Nguess(kds, kdr, pisr, pirs, Nin, kobs)
dt = 0.001;
t_steps = 7./dt;
t_int = 7;% interval is 7 days

 
% Sets up transition matrix

% diagonals must add to kobs
% assume that sensitive cells are dying while resistant cells are able to
% proliferate
% we have the equation e^kdst + e^kdrt = e^kt

% Unconstrained k values ( arbitrary and not related to k
k_init_guess = zeros(2,2); 
k_init_guess(1,1) = kds; % LHS diagonal kds
%k_init_guess(2,2) = (sum(Nin)./Nin(2))*(k_wk-k_init_guess(1,1)) + k_init_guess(1,1); % kdr scaled properly
k_init_guess(2,2) = kdr;
% Note that there are only three true unknowns here: pi s->r,pi r->s, and kds
% Therefore these are the paramters we need to optimize.



pi_init_guess = 0.5*ones(2,2);
pi_init_guess(1,2) = pisr; % pi s->r
pi_init_guess(1,1) = 1-pi_init_guess(1,2); % pi s->s ****
pi_init_guess(2,1) = pirs; % pi r->s
pi_init_guess(2,2) = 1-pi_init_guess(2,1); % pi r->r

% diagonals must add to kobs
% assume that sensitive cells are dying while resistant cells are able to
% proliferate
% RETURN TO THIS IF FDM DOES NOT WORK
% k_wk = kobs(1).^7; % observed k over a week
% k_init_guess = zeros(2,2); 
% k_init_guess(1,1) = kds; % LHS diagonal kds
% k_init_guess(2,2) = (sum(Nin)./Nin(2))*(k_wk-k_init_guess(1,1)) + k_init_guess(1,1); % kdr scaled properly

% Note that there are now four unknowns here: pi s->r,pi r->s, kds, and kdr
% Therefore these are the paramters we need to optimize.

% Now find the first estimate of Nout
Nout_guess_prev = Nin*(k_init_guess*pi_init_guess);

N_guess = zeros(t_steps, 2);
N_guess(1,:) = Nin;



% Allow for state transitions
for t = 2:t_steps
    N_guess_curr = N_guess(t-1,:)';
    N_guess(t,:) = N_guess_curr + dt*((k_init_guess*pi_init_guess)*(N_guess_curr));
    N_guess(t,:) = N_guess(t,:)';
end

Nout_guess_fin = N_guess(t_steps, :); % outputs final guess 


end