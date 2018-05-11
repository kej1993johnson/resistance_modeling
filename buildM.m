function [M] = buildM(sx)

% sx = # of points in the x-direction




A = zeros(sx,sx);
v = zeros(1, sx);
r = ones(1,sx);
A = diag(v);
A = diag(v,0);
b = diag(r,1);
c = diag(r,-1);


M = A + b(1:sx, 1:sx) + c(1:sx, 1:sx);

%A = sparse(A);
end