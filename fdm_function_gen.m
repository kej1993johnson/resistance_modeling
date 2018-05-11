function [ N] = fdm_function_gen( t, sx, sy, dt, dx, k, dy, N, D, carcap, dD_dx, dD_dy);
%FDM Function: This function takes the input of the current N function and
%D function an updates the new N at time t+ dt using central difference
%approximation in spaces and forward difference approximation in time.
%This function is to be called for all x and y not equal to the boundary
%conditions
% D = 0.5.*D0((x./100).^2 + (y./100).^2);
% dD_dx = (D0.*x)./100;
% dD_dy = (D0.*y)./100;
% iterate through all xs at each y (could also be all ys at each x, order
% shouldn't matter. This fills the entire 2D N (x, y) for a single time
% point that is inputted to the function
for x = 2: (sx -1)
for y = 2: (sy -1)
   first_deriv_x = (1./(2.*dx).*dD_dx(x,y)).*(N(x+1,y,t) - N(x-1,y,t));
   second_deriv_x = (D(x,y)./(dx.^2)).*(N(x+1,y,t)-2.*N(x,y,t) + N(x-1,y,t));
   first_deriv_y = (1./(2.*dy).*dD_dy(x,y)).*(N(x,y+1,t) - N(x,y-1,t));
   second_deriv_y = (D(x,y)./(dy.^2)).*(N(x,y+1,t)-2.*N(x,y,t) + N(x,y-1,t));
   k_term = k.*N(x,y,t).*(1- (N(x,y,t)./carcap));
   
   N(x,y,t+1) = N(x,y,t) + dt.*( first_deriv_x + second_deriv_x + first_deriv_y + second_deriv_y + k_term);
    
%     N(x,y,t+1) = N(x,y,t)... 
% + dt.*(((1./(2.*dx)).*((D0.*x.*dx)./100).*(N(x+1,y,t) - N(x-1,y,t))...
% +(1./(2.*dy)).*((D0.*y.*dy)./100).* (N(x,y+1,t) - N(x,y-1,t))...
% + (0.5.*D0.*((x.*dx./100).^2 + (y.*dy./100).^2)./(dx.^2))...
% .*(N(x+1, y, t) + N(x-1,y,t))...
% + (0.5.*D0.*((x.*dx./100).^2 + (y.*dy./100).^2)./(dy.^2))...
% .*(N(x,y+1,t)+N(x,y-1,t))...
% + N(x,y,t).*k.*(1- (N(x,y,t)./carcap))));

end
end

end

