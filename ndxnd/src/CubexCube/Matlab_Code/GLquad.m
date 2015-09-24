function [x,w]=GLquad(n,a,b)

% computes nodes x and weights w of Gauss-Legendre Quadrature Rule on [a,b]
% by solving a symmetric tridiagonal eigenvalue problem (Golub/Welsch 1969)

if nargin==1 
    a=0; b=1; 
end

[x,w]=gauleg(n);

x=x'*(b-a)/2+(b+a)/2;
w=w'*(b-a)/2;

end

