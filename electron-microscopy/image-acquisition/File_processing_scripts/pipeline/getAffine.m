function [A] = getAffine(x, y, xp, yp)
% Calculates the affine transform based on the corresponding points
% estimates (x,y,1) -> (xp, yp,1)
% Returns 3x3 matrix that represents the affine transform so that
% A*(x,y,1)' = (xp, yp, 1);
%
% Aleksandar Zlateski, 2009

[m n] = size(x);
A = zeros(2*n,6);
b = zeros(2*n,1);

for i=1:n,
    A(2*i-1,:)  = [x(i),y(i),1,0,0,0];
    A(2*i,:)    = [0,0,0,x(i),y(i),1];
    b(2*i-1)    = xp(i);
    b(2*i)      = yp(i);
end;

H = A\b;
A = [H(1), H(2), H(3);
     H(4), H(5), H(6);
     0   , 0   , 1];
