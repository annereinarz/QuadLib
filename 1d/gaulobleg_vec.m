% GAULOB  computes Gauss-Lobatto quadrature points in [-1,1]
%
% function [x,w]=gaulob(n)
%
%input: n = Anzahl Quadraturpunkte
%ouput: x = Vektor mit den Stuetzstellen, 
%       w = Vektor mit den zugehoerigen Gewichten
%function [x,w,count]=gaulobleg_vec(n)
function [x,w]=gaulobleg_vec(n)

m=floor((n-1)/2);
z=cos(pi*(1:m)/(n-1));
dz=1;
count = 0;
while max(abs(dz))>eps    
    b=1;    
    c=z; 
    for k=2:n-1
        a=b;
        b=c;
        c=((2*k-1)*z.*b-(k-1)*a)/k;
    end
    %dz=(z.*c-b)./(n*c+2*z.*(z.*c-b)./(1-z.^2));
    dz=(z.*c-b)./(n*c);
    z=z-dz; 
    count = count+1;
end
x=[1 z]';
x(n:-1:n-m) = -x(1:m+1);
w=2./((n-1)*n*[1 c.^2]');
w(n:-1:n-m) = w(1:m+1);
if n>2
   x=x';
end
end
