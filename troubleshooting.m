% troubleshooting eqn
x = 46
y = 46
t = 1
N(x,y,t+1) = N(x,y,t)... 
+ dt.*(((1./(2.*dx)).*((D0.*x.*dx)./100).*(N(x+1,y,t) - N(x-1,y,t))...
+(1./(2.*dy)).*((D0.*y.*dy)./100).* (N(x,y+1,t) - N(x,y-1,t))...
+ (0.5.*D0.*((x.*dx./100).^2 + (y.*dy./100).^2)./(dx.^2))...
.*(N(x+1, y, t) + N(x-1,y,t))...
+ (0.5.*D0.*((x.*dx./100).^2 + (y.*dy./100).^2)./(dy.^2))...
.*(N(x,y+1,t)+N(x,y-1,t))...
+ N(x,y,t).*k.*(1- (N(x,y,t)./carcap))));
N(46,46, 2)
N(46,46,1)

%%
% for x = 2: (sx -1)
% for y = 2: (sy -1)
x = 45
y = 45
t = 1
dD_dx(5,10)
for t = 1:10
for x = 45:55
    for y = 45:55
   first_deriv_x = (1./(2.*dx).*dD_dx(x,y)).*(N(x+1,y,t) - N(x-1,y,t));
   second_deriv_x = (D(x,y)./(dx.^2)).*(N(x+1,y,t)-2.*N(x,y,t) + N(x-1,y,t));
   first_deriv_y = (1./(2.*dy).*dD_dy(x,y)).*(N(x,y+1,t) - N(x,y-1,t));
   second_deriv_y = (D(x,y)./(dy.^2)).*(N(x,y+1,t)-2.*N(x,y,t) + N(x,y-1,t));
   k_term = k.*N(x,y,t).*(1- (N(x,y,t)./carcap));
   
   N(x,y,t+1) = N(x,y,t) + dt.*( first_deriv_x + second_deriv_x + first_deriv_y + second_deriv_y + k_term);
    end
end
end
N(47,47,6)
   N(45, 45, 4) 
%     N(x,y,t+1) = N(x,y,t)... 
% + dt.*(((1./(2.*dx)).*((D0.*x.*dx)./100).*(N(x+1,y,t) - N(x-1,y,t))...
% +(1./(2.*dy)).*((D0.*y.*dy)./100).* (N(x,y+1,t) - N(x,y-1,t))...
% + (0.5.*D0.*((x.*dx./100).^2 + (y.*dy./100).^2)./(dx.^2))...
% .*(N(x+1, y, t) + N(x-1,y,t))...
% + (0.5.*D0.*((x.*dx./100).^2 + (y.*dy./100).^2)./(dy.^2))...
% .*(N(x,y+1,t)+N(x,y-1,t))...
% + N(x,y,t).*k.*(1- (N(x,y,t)./carcap))));

% end
% end
% 
% end

