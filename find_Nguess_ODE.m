function [Nout_guess_fin] = find_Nguess_ODE(gs,gr,ksr,krs, Nin, kmark)
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


end