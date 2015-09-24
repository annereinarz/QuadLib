function [x,w]=GLquad(n,a,b)

% function [x,w]=GLquad(n,a,b) 
%
% computes nodes x and weights w of Gauss-Legendre Quadrature Rule on [a,b]
% by solving a symmetric tridiagonal eigenwalue problem (Golub/Welsch 1969)
% uses GLquad_sym
% 
% x = quardature nodes
% w = quadrature weights
% n = number of quadrature pairs
%
% (C) Alexey Chernov, ETH Zurich, 2008

if nargin==1 
    a=0; b=1; 
end

%[x,w]=GLquad_sym(n);
[x,w]=gauleg(n);

x=x'*(b-a)/2+(b+a)/2;
w=w'*(b-a)/2;

end

function [x,w]=GLquad_sym(n)

% function [x,w]=GLquad_sym(n) 
%
% computes nodes x and weights w of Gauss-Legendre Quadrature Rule on [-1,1]
% by solving a symmetric tridiagonal eigenwalue problem (Golub/Welsch 1969)
% 
% x = quardature nodes
% w = quadrature weights
% n = number of quadrature pairs
%
% (C) Alexey Chernov, ETH Zurich, 2008

b=1:n-1; 
b=b./sqrt(4*b.*b-1);
J=diag(b,-1)+diag(b,1); 
[evec,eval]=eig(J);

% normalize eigenvectors

for i=1:n
    evec(:,i)=evec(:,i)./norm(evec(:,i));
end

x=diag(eval);
w=(2*(evec(1,:).*evec(1,:)))';

end

function [xg,wg] = gauleg(n)
%      [xg,wg] = gauleg(n)
%      find Gauss nodes xg(1),...xg(n) and weights wg(1),...,wg(n) for [-1,1]
%
%      Then integral of f(x) over x in [-1,1]
%      can be approximated by   wg(1)*f(xg(1)) + ... + wg(n)*f(xg(n))
%
eps = 1.5*2^(-53);
m=(n+1)/2;
for i=1:m
  z=cos(pi*(i-0.25)/(n+0.5));
  while 1
    p1=1;
    p2=0;
    for j=1:n
      p3=p2;
      p2=p1;
      p1=((2*j-1)*z*p2-(j-1)*p3)/j;
    end
    pp=n*(z*p1-p2)/(z*z-1);
    z1=z;
    z=z1-p1/pp;
    if abs(z-z1) <= eps
      break
    end
  end
  xg(i) = -z;
  xg(n+1-i) = z;
  wg(i) = 2/((1-z*z)*pp*pp);
  wg(n+1-i)=wg(i);
end

end
