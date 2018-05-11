close all, clear all, clc
% This script solves the equation du/dt = -du/dx in 4 ways:
% 1. Leap-Frog Method
% 2. Lax-Friendrichs Method
% 3. Upwinding Method
% 4. MUSCL Method



dt = 0.01;
T1= 1;
T2 = 5;
T3 = 10; % first-time point
h = 0.25; % first mesh spacing
domain_x = -5:h:15; % 1x 20
sx = length(domain_x); % length of x domain in grid points

u1 = zeros(sx,1); %
u2 = zeros(sx,1);

%% First Initial Condition
% First find u(:,1) and u(:,2) for first two time points using the true
% solution
% First initial condition of u(x,0) = x
for x = 1:sx
    % Since domain is shifted, need to adjust for this, so true first x =
    % -5
        u1(x,1) = (x-1)*h-5;
        u1(x,2) =((u1(x,1)-dt));
end
%% Second Initial Condition
for x = 1:((sx-1)*0.25) % waant to run for the first 5
    % Since domain is shifted, need to adjust for this, so true first x =
    % -5
        u2(x,1) = 1;
end
for x =((sx-1)*0.25)+1:sx
    % Since domain is shifted, need to adjust for this, so true first x =
    % -5
        u2(x,1) = 0;
end

for x =1:sx
    u2(x,2) = (((x-1)*h)-5)-dt;
end
%% 1. Leap-Frog Method
% Have already created domains for u where each column is an x point and
% each row is a time point of dt
Q = buildQ(sx, dt, h);

for t = 2:((T3/dt)-1)
    u1(1,t+1) = ((-5-((t)*dt))); % set boundary at -5 = true solution (-5-t))
    u1(sx,t+1) = ((15-((t)*dt))); % set boundary at 15 = true solution (15-t))
    u1(2:sx-1,t+1) = u1(2:sx-1,t-1) + Q(2:sx-1, 1:sx)*u1(1:sx,t); % only find inner x values
    u2(1,t+1) = ((-5-((t)*dt))); % set boundary at -5 = true solution (-5-t)
    u2(sx,t+1) = ((15-((t)*dt)));% set boundary at 15 = true solution (15-t)
    u2(2:sx-1,t+1) = u2(2:sx-1,t-1) + Q(2:sx-1, 1:sx)*u2(1:sx,t);
    
    
end
%%
figure(1)
subplot(2,3,1)
set(gca,'FontSize',20,'LineWidth',3)
plot(domain_x, u1(:,T1/dt))
xlim([-5,15])
xlabel('x')
ylabel('u')
title('Solution T=1s, u0=x')
subplot(2,3,2)
set(gca,'FontSize',20,'LineWidth',3)
plot(domain_x, u1(:,T2/dt))
xlabel('x')
xlim([-5,15])
ylabel('u')
title('Solution T=5s, u0=x')
subplot(2,3,3)
set(gca,'FontSize',20,'LineWidth',3)
plot(domain_x, u1(:,T3/dt))
xlabel('x')
xlim([-5,15])
ylabel('u')
title('Solution T=10s, u0=x')

subplot(2,3,4)
set(gca,'FontSize',20,'LineWidth',3)
plot(domain_x, u2(:,T1/dt))
xlim([-5,15])
xlabel('x')
ylabel('u')
title('Solution T=1s, u0=0,1')
subplot(2,3,5)
set(gca,'FontSize',20,'LineWidth',3)
plot(domain_x, u2(:,T2/dt))
xlabel('x')
xlim([-5,15])
ylabel('u')
title('Solution T=5s, u0=0,1')
subplot(2,3,6)
set(gca,'FontSize',20,'LineWidth',3)
plot(domain_x, u2(:,T3/dt))
xlabel('x')
xlim([-5,15])
ylabel('u')
title('Solution T=10s, u0=0,1')

% This is the solution for a mesh of 80 points and a dt = 0.1. Appears to
% become unstable at after 1s.... not sure if this is due to coding error
% in boundary conditions/implementation or if this is instability of the
% method...

%% Lax-Friendrich's Method
% reset  everything beyond first time point  in u1 and u2 to zero
u1 = u1(:,1);
u2 = u2(:,1);

Q = buildQ(sx, dt, h);
M = buildM(sx);

for t = 1:((T3/dt)-1)
    u1(1,t+1) = ((-5-((t)*dt))); % set boundary at -5 = true solution (-5-t))
    u1(sx,t+1) = ((15-((t)*dt))); % set boundary at 15 = true solution (15-t))
    u1(2:sx-1,t+1) = 0.5.*M(2:sx-1, 1:sx)*u1(1:sx,t)- 0.5.*(Q(2:sx-1, 1:sx)*u1(1:sx,t)); % only find inner x values
    u2(1,t+1) = ((-5-((t)*dt))); % set boundary at -5 = true solution (-5-t))
    u2(sx,t+1) = ((15-((t)*dt))); % set boundary at 15 = true solution (15-t))
    u2(2:sx-1,t+1) = 0.5.*(M(2:sx-1, 1:sx)*u2(1:sx,t))- 0.5.*(Q(2:sx-1, 1:sx)*u2(1:sx,t)); % only find inner x values

    
end
%% Upwinding method
% reset initial conditions
u1 = u1(:,1);
u2 = u2(:,1);

Q_u = buildQ_u(sx,dt,h);


for t = 1:((T3/dt)-1)
    u1(1,t+1) = ((-5-((t)*dt))); % set boundary at -5 = true solution (-5-t))
    u1(sx,t+1) = ((15-((t)*dt))); % set boundary at 15 = true solution (15-t))
    u1(2:sx-1,t+1) = Q_u(2:sx-1, 1:sx)*u1(1:sx,t); % only find inner x values
    u2(1,t+1) = ((-5-((t)*dt))); % set boundary at -5 = true solution (-5-t))
    u2(sx,t+1) = ((15-((t)*dt))); % set boundary at 15 = true solution (15-t))
    u2(2:sx-1,t+1) =Q_u(2:sx-1, 1:sx)*u2(1:sx,t); % only find inner x values

    
end
%% MUSCL method
% will have to write a code to calculate f(Unj) depending on the difference
% between intervals
u1 = u1(:,1);
u2 = u2(:,1);


for t = 1:((T3/dt)-1)
    u1(1,t+1) = ((-5-((t)*dt))); % set boundary at -5 = true solution (-5-t))
    u1(2, t+1) = ((-(5-h)-((t)*dt))); % set boundary at -4.9 = true solution (-4.9-t)
    u1(sx,t+1) = ((15-((t)*dt))); % set boundary at 15 = true solution (15-t))
    u1(sx-1,t+1) = (((15-h)-((t)*dt)));% set boundary at 14.9 = true solution (14.9-t)
    for x = 3:sx-3
        
        ULplus = -1/6*u1(x-1,t)+5/6*u1(x,t)-1/3*u1(x+1,t);
        URplus = 1/3*u1(x,t)+5/6*u1(x+1,t)-1/6*u1(x+2,t);
        ULminus = -1/6*u1(x-2,t)+5/6*u1(x-1,t)-1/3*u1(x,t);
        URminus = 1/3*u1(x-1,t)+5/6*u1(x,t)-1/6*u1(x+1,t);
        % Slope limiter calculation
        if ULplus <URplus
            F_plus = ULplus;
        end
        if URplus< ULplus
            F_plus = URplus;
        end
        
         if ULminus <URminus
            F_minus = ULminus;
        end
        if URminus< ULminus
            F_minus = URminus;
        end
        if URplus*ULplus <= 0
            F_plus = 0;
        end
        if URminus*ULminus <= 0
            F_minus = 0;
        end
        
        u1(x,t+1) = u1(x, t) - dt/h*(F_plus-F_minus);
    end
    
end
%%
for t = 1:((T3/dt)-1)

    u2(1,t+1) = ((-5-((t)*dt))); % set boundary at -5 = true solution (-5-t)
     u2(2, t+1) = ((-(5-h)-((t)*dt))); % set boundary at -4.9 = true solution (-4.9-t)
    u2(sx,t+1) = ((15-((t)*dt)));% set boundary at 15 = true solution (15-t)
     u2(sx-1,t+1) = (((15-h)-((t)*dt))); % % set boundary at 14.9 = true solution (14.9-t)
    for x = 3:sx-3
        
        ULplus = -1/6*u2(x-1,t)+5/6*u2(x,t)-1/3*u2(x+1,t);
        URplus = 1/3*u2(x,t)+5/6*u2(x+1,t)-1/6*u2(x+2,t);
        ULminus = -1/6*u2(x-2,t)+5/6*u2(x-1,t)-1/3*u2(x,t);
        URminus = 1/3*u2(x-1,t)+5/6*u2(x,t)-1/6*u2(x+1,t);
        % Slope limiter calculation
        if ULplus <URplus
            F_plus = ULplus;
        end
        if URplus< ULplus
            F_plus = URplus;
        end
        
         if ULminus <URminus
            F_minus = ULminus;
        end
        if URminus< ULminus
            F_minus = URminus;
        end
        if URplus*ULplus <= 0
            F_plus = 0;
        end
        if URminus*ULminus <= 0
            F_minus = 0;
        end
        
        u2(x,t+1) = u2(x, t) - dt/h*(F_plus-F_minus);
    end
    
    
end






