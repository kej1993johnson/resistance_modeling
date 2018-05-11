function [ N(:,:,t+1) ] = fdm_function_gen(t, sx, sy, dx, dy, D0, N(:,:,t), D(:,:));
%FDM Function: This function takes the input of the current N function and
%D function an updates the new N at time t+ dt using central difference
%approximation in spaces and forward difference approximation in time.
%This function is to be called for all x and y not equal to the boundary
%conditions
D = 0.5.*D0((x./100).^2 + (y./100).^2);
dD_dx = (D0.*x)./100;
dD_dy = (D0.*y)./100;
% iterate through all xs at each y (could also be all ys at each x, order
% shouldn't matter. This fills the entire 2D N (x, y) for a single time
% point that is inputted to the function
for x = .1: ((sx./dx) -dx)
for y = .1: ((sy./dy) -dy)
N(x,y,t+1) = N(x,y,t) + dt((1./(2.*dx)).*((D0.*x)./100).*(N(x+1,y,t) - N(x-1,y,t))+ (1./(2.*dy)).*((D0.*y)./100).* (N(x,y+1,t) - N(x,y-1,t))+ (D./(dx.^2)).*(N(x+1, y, t) + N(x-1,y,t)) + (0.5.*D0((x./100).^2 + (y./100).^2)./(dy.^2)).*(N(x, y+1, t)+N(x, y-1, t)) + N(x, y, t)).* k.*(1- (N(x, y, t)./carcap));
end
end

end

