% Johnson_K_hw_1.m
% Created on: February 2nd, 2017
% Created by: Kaitlyn Johnson
%Please include comments throughout your work for partial credit.


%% Problem 3a
clear all; clc; 
% Initialize your parameters and simulation domain
D0 = 0.05;
k = 0; % for 3a
carcap = 100;
dt = 0.01;
num_iters = 500; % note this is days observed./dt
% Want to evaluate in 10 x 10 at .1 spacial step
dx = 0.1;
dy = 0.1;
sx = 100; % number of spatial grids in the x direction (10/ .1)
sy = 100; % number of spatial grids in the y direction (10/ .1)

% need to figure out how to get 3D geometry?
%% Define Domain
% N(x,y,t)-> N is a function of x, y ,t where x and y range from x = 0 to x
% =10 and y = 0 to y = 10. Start by defining this domain as all zeros, then
% initialize with 101 x 101 grids and 500 time steps (each for .01 days over 5
% days). Chose to use 101 x 101 so that when performing center difference
% N(x-1) and N(y-1) are not calling indices = 0
N=zeros([102,102,500]);
% D(x, y) -> D is a function of x and y ranging from x = 0 to x = 10.
% Initialize the diffusion coefficient at each spatial point
D= zeros([101,101]);
% Make diffusion matrix since it is constant throughout time. i and j terms
% have minus one because of shift imposed by padding the field of view.
% Multiply i and j by x and y because each grid is .1 mm and we assume the
% diffusion equation is given in terms of millimeters
for i = 2:101
    for j = 2:101
        D(i, j) = (D0./2).*((((i-1)./100).^2) + (((j-1)./100).^2));
    end    
end
% This sets D
% Now set D_dx = the derivative of D taken with respect to x
for i = 2:101
    for j = 2:101
        dD_dx(i, j) = D0.*(i-1)./(100.^2);
        dD_dy(i, j) = D0.*(j-1)./(100.^2);
    end    
end
% Now let these be inputs to function
% Images to check
imagesc(dD_dx(:,:))
max(max(D));
%%
% Now initialize beginnign of tumor growth with a 1 mm (10 grid edge)
N(45:55, 45:55, 1) = 0.75.* carcap;


for t = 1:(num_iters)
    for x = 2:101
    for y = 2:101
   first_deriv_x = (1./(2.*dx).*dD_dx(x,y)).*(N(x+1,y,t) - N(x-1,y,t));
   second_deriv_x = (D(x,y)./(dx.^2)).*(N(x+1,y,t)-2.*N(x,y,t) + N(x-1,y,t));
   first_deriv_y = (1./(2.*dy).*dD_dy(x,y)).*(N(x,y+1,t) - N(x,y-1,t));
   second_deriv_y = (D(x,y)./(dy.^2)).*(N(x,y+1,t)-2.*N(x,y,t) + N(x,y-1,t));
   k_term = k.*N(x,y,t).*(1- (N(x,y,t)./carcap));
   N(x,y,t+1) = N(x,y,t) + dt.*( first_deriv_x + second_deriv_x + first_deriv_y + second_deriv_y + k_term);
   
   % Now apply boundary conditions to edges. Since we want the derivatives
   % to be 0 at the edge, set the left and right terms of the edge to be
   % equal so that the first_derivatives go to 0 at the edges
   N(1,y,t) = N(3, y, t); % makes first_deriv_x = 0 at x = 2 which is left boundary
   N(100, y, t+1) = N(102, y, t); % makes first_deriv_x =0 at x = 101 which is right boundary
   N(x,1,t) = N(x, 3, t); % makes first_deriv_y = 0 at y = 2 which is left boundary
   N(x, 100, t) = N(x, 102, t); % makes first_deriv_x = 0 at x = 101 which is right boundary 
   
   % Now apply boundary conditions to the corners. At the corners both
   % first derivatives will be equal to 0, and both second derivatives will
   % be equal to 0 (since the derivative of 0 is 0). Therefore the only
   % term present is the k term at the corners. We can write this as
   N(2,2,t+1) = k_term;
   N(2,101, t+1) = k_term;
   N(101,101,t+1) = k_term;
   N(101, 2, t+1) = k_term;
  
  
 end
end
Ntotal(t,1) = sum(sum(N(:,:,t)));
end



%%
imagesc(N(:,:,500))

%%
% Problem 3a plot
Nchange(:,1) = Ntotal(:,1) - Ntotal(1,1);
% Display Results
figure(1)
plot(1:500, Nchange(:,1),'LineWidth',3)
xlim([1,400])
xlabel('Iteration','FontSize',20)
ylabel('Change in Cell Number','FontSize',20)
title('Problem 3a, k =0, t=5 days, 500 iterations', 'FontSize',14)
set(gca,'LineWidth',1.5,'FontSize',20)

% hold off
% imagesc(N(:,:,1));
% imagesc(N(:,:,500));

%% Problem 3b
clear all; clc; 
% Initialize your parameters and simulation domain
D0 = 0.05;
k = 0; % for 3a
carcap = 100;
dt = 0.01;
num_iters = 400; % note this is days observed./dt
% Want to evaluate in 10 x 10 at .1 spacial step
dx = 0.1;
dy = 0.1;

% Define Domain
N=zeros([102,102,400]); % Now only 400 iterations
D= zeros([101,101]);

for i = 2:101
    for j = 2:101
        D(i, j) = (D0./2).*((((i-1)./100).^2) + (((j-1)./100).^2));
    end    
end
% This sets D
for i = 2:101
    for j = 2:101
        dD_dx(i, j) = D0.*(i-1)./(100.^2);
        dD_dy(i, j) = D0.*(j-1)./(100.^2);
    end    
end

% Now initialize beginning of tumor growth with a 1 mm (10 grid edge)
N(45:55, 45:55, 1) = 0.75.* carcap;

for m = 1:6
    dt = .01.*(2.*m-1); % adds another loop that runs through .02 increases in dt
    
for t = 1:(num_iters)
    for x = 2:101
    for y = 2:101
   first_deriv_x = (1./(2.*dx).*dD_dx(x,y)).*(N(x+1,y,t) - N(x-1,y,t));
   second_deriv_x = (D(x,y)./(dx.^2)).*(N(x+1,y,t)-2.*N(x,y,t) + N(x-1,y,t));
   first_deriv_y = (1./(2.*dy).*dD_dy(x,y)).*(N(x,y+1,t) - N(x,y-1,t));
   second_deriv_y = (D(x,y)./(dy.^2)).*(N(x,y+1,t)-2.*N(x,y,t) + N(x,y-1,t));
   k_term = k.*N(x,y,t).*(1- (N(x,y,t)./carcap));
   N(x,y,t+1) = N(x,y,t) + dt.*( first_deriv_x + second_deriv_x + first_deriv_y + second_deriv_y + k_term);
   
   N(1,y,t) = N(3, y, t); % makes first_deriv_x = 0 at x = 2 which is left boundary
   N(100, y, t+1) = N(102, y, t); % makes first_deriv_x =0 at x = 101 which is right boundary
   N(x,1,t) = N(x, 3, t); % makes first_deriv_y = 0 at y = 2 which is left boundary
   N(x, 100, t) = N(x, 102, t); % makes first_deriv_x = 0 at x = 101 which is right boundary 
   
   N(2,2,t+1) = k_term;
   N(2,101, t +1) = k_term;
   N(101,101,t+1) = k_term;
   N(101, 2, t+1) = k_term;
  
 end
    end
Ntotal(t,1) = sum(sum(N(:,:,t))); % at end of time calculates Ntotal
end

Nchange(:,m) = Ntotal(:,1) - Ntotal(1,1); % put change in cell number for each iteration 
%                                           % with each column containing a
%                                           % diffferent delta t vale

end


% Problem 3b plot

% Display Results
figure(2)
hold off
%plot(1:400, Nchange(:,1),'LineWidth',3)
% % hold on
plot(1:400, Nchange(:,6),'LineWidth',3)
% plot(1:400, Nchange(:,6),'LineWidth',3)
xlim([1,400])
xlabel('Iteration','FontSize',20)
ylabel('Change in Cell Number','FontSize',20)
title('Problem 3b, k =0, dt = 0.11', 'FontSize',14)
set(gca,'LineWidth',1.5,'FontSize',20)

%% Problem 3c
clear all; clc; 

% Initialize your parameters and simulation domain
for m = 1:7
D0 = 0.0015.*(2.^(m-1));
k = 0.003; % for 3a
carcap = 100;
dt = 0.01;
num_iters = 400; % note this is days observed./dt
% Want to evaluate in 10 x 10 at .1 spacial step
dx = 0.1;
dy = 0.1;

% Define Domain
N=zeros([102,102,400]); % Now only 400 iterations
D= zeros([101,101]);

for i = 2:101
    for j = 2:101
        D(i, j) = (D0./2).*((((i-1)./100).^2) + (((j-1)./100).^2));
    end    
end
% This sets D
for i = 2:101
    for j = 2:101
        dD_dx(i, j) = D0.*(i-1)./(100.^2);
        dD_dy(i, j) = D0.*(j-1)./(100.^2);
    end    
end

% Now initialize beginning of tumor growth with a 1 mm (10 grid edge)
N(45:55, 45:55, 1) = 0.75.* carcap;

for t = 1:(num_iters)
    for x = 2:101
    for y = 2:101
   first_deriv_x = (1./(2.*dx).*dD_dx(x,y)).*(N(x+1,y,t) - N(x-1,y,t));
   second_deriv_x = (D(x,y)./(dx.^2)).*(N(x+1,y,t)-2.*N(x,y,t) + N(x-1,y,t));
   first_deriv_y = (1./(2.*dy).*dD_dy(x,y)).*(N(x,y+1,t) - N(x,y-1,t));
   second_deriv_y = (D(x,y)./(dy.^2)).*(N(x,y+1,t)-2.*N(x,y,t) + N(x,y-1,t));
   k_term = k.*N(x,y,t).*(1- (N(x,y,t)./carcap));
   N(x,y,t+1) = N(x,y,t) + dt.*( first_deriv_x + second_deriv_x + first_deriv_y + second_deriv_y + k_term);
   
   N(1,y,t) = N(3, y, t); % makes first_deriv_x = 0 at x = 2 which is left boundary
   N(100, y, t+1) = N(102, y, t); % makes first_deriv_x =0 at x = 101 which is right boundary
   N(x,1,t) = N(x, 3, t); % makes first_deriv_y = 0 at y = 2 which is left boundary
   N(x, 100, t) = N(x, 102, t); % makes first_deriv_x = 0 at x = 101 which is right boundary 
   
   N(2,2,t+1) = k_term;
   N(2,101, t +1) = k_term;
   N(101,101,t+1) = k_term;
   N(101, 2, t+1) = k_term;
  
 end
    end
Ntotal(t,1) = sum(sum(N(:,:,t))); % at end of time calculates Ntotal
end

Nmid(:, m) = N(:, 51, 400);
end

% Problem 3b plot
%%
% Display Results
figure(3)
hold off
plot(1:102, Nmid(:,1), 'LineWidth', 2)
hold on
plot(1:102, Nmid(:,2), 'LineWidth', 2)
hold on
plot(1:102, Nmid(:,3), 'LineWidth', 2)
plot(1:102, Nmid(:,4), 'LineWidth', 2)
plot(1:102, Nmid(:,5), 'LineWidth', 2)
plot(1:102, Nmid(:,6), 'LineWidth', 2)
plot(1:102, Nmid(:,7),'LineWidth', 2)
xlabel('Tumor Dimension (mm)','FontSize',20)
ylabel('Number of Tumor Cells','FontSize',20)
title('Problem 3c, k =0.003, increasing diffusion', 'FontSize',14)
legend('D/k = 0.5', 'D/k = 1', 'D/k = 2','D/k = 4', 'D/k = 8', 'D/k = 16', 'D/k = 32')
set(gca,'LineWidth',1.5,'FontSize',20)

%% Problem 3d
clear all; clc; 
clear all; clc; 
% Initialize your parameters and simulation domain
D0 = 0.05;
k = 2.5; % for 3a
carcap = 100;
dt = 0.01;
num_iters = 500; % note this is days observed./dt
% Want to evaluate in 10 x 10 at .1 spacial step
dx = 0.1;
dy = 0.1;
sx = 100; % number of spatial grids in the x direction (10/ .1)
sy = 100; % number of spatial grids in the y direction (10/ .1)

% need to figure out how to get 3D geometry?
% Define Domain
% N(x,y,t)-> N is a function of x, y ,t where x and y range from x = 0 to x
% =10 and y = 0 to y = 10. Start by defining this domain as all zeros, then
% initialize with 101 x 101 grids and 500 time steps (each for .01 days over 5
% days). Chose to use 101 x 101 so that when performing center difference
% N(x-1) and N(y-1) are not calling indices = 0
N=zeros([102,102,500]);
% D(x, y) -> D is a function of x and y ranging from x = 0 to x = 10.
% Initialize the diffusion coefficient at each spatial point
D= zeros([101,101]);
% Make diffusion matrix since it is constant throughout time. i and j terms
% have minus one because of shift imposed by padding the field of view.
% Multiply i and j by x and y because each grid is .1 mm and we assume the
% diffusion equation is given in terms of millimeters
for i = 2:101
    for j = 2:101
        D(i, j) = (D0./2).*((((i-1)./100).^2) + (((j-1)./100).^2));
    end    
end
% This sets D
% Now set D_dx = the derivative of D taken with respect to x
for i = 2:101
    for j = 2:101
        dD_dx(i, j) = D0.*(i-1)./(100.^2);
        dD_dy(i, j) = D0.*(j-1)./(100.^2);
    end    
end
% Now let these be inputs to function
% Images to check
imagesc(dD_dx(:,:))
max(max(D));

% Now initialize beginnign of tumor growth with a 1 mm (10 grid edge)
N(45:55, 45:55, 1) = 0.75.* carcap;


for t = 1:(num_iters)
    for x = 2:101
    for y = 2:101
   first_deriv_x = (1./(2.*dx).*dD_dx(x,y)).*(N(x+1,y,t) - N(x-1,y,t));
   second_deriv_x = (D(x,y)./(dx.^2)).*(N(x+1,y,t)-2.*N(x,y,t) + N(x-1,y,t));
   first_deriv_y = (1./(2.*dy).*dD_dy(x,y)).*(N(x,y+1,t) - N(x,y-1,t));
   second_deriv_y = (D(x,y)./(dy.^2)).*(N(x,y+1,t)-2.*N(x,y,t) + N(x,y-1,t));
   k_term = k.*N(x,y,t).*(1- (N(x,y,t)./carcap));
   N(x,y,t+1) = N(x,y,t) + dt.*( first_deriv_x + second_deriv_x + first_deriv_y + second_deriv_y + k_term);
   
   % Now apply boundary conditions to edges. Since we want the derivatives
   % to be 0 at the edge, set the left and right terms of the edge to be
   % equal so that the first_derivatives go to 0 at the edges
   N(1,y,t) = N(3, y, t); % makes first_deriv_x = 0 at x = 2 which is left boundary
   N(100, y, t+1) = N(102, y, t); % makes first_deriv_x =0 at x = 101 which is right boundary
   N(x, 1,t) = N(x, 3, t); % makes first_deriv_y = 0 at y = 2 which is left boundary
   N(x, 100, t) = N(x, 102, t); % makes first_deriv_x = 0 at x = 101 which is right boundary 
   
   % Now apply boundary conditions to the corners. At the corners both
   % first derivatives will be equal to 0, and both second derivatives will
   % be equal to 0 (since the derivative of 0 is 0). Therefore the only
   % term present is the k term at the corners. We can write this as
   N(2,2,t+1) = k_term;
   N(2,101, t+1) = k_term;
   N(101,101,t+1) = k_term;
   N(101, 2, t+1) = k_term;
  
  
 end
end
Ntotal(t,1) = sum(sum(N(:,:,t)));
end
%%
figure (4)
hold off
imagesc(N(:,:,100))
xlabel('X dimension (0.1 mm)','FontSize',20)
ylabel('Y-dimension(0.1 mm)','FontSize',20)
title('Problem 3d, k =2.5, t=1 day, 500 iterations', 'FontSize',14)
set(gca,'LineWidth',1.5,'FontSize',20)
