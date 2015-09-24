% GAULOB  computes Gauss-Lobatto quadrature points in [-1,1]
%
% function [x,w]=gaulob(n)
%
%input: n = Anzahl Quadraturpunkte
%ouput: x = Vektor mit den Stuetzstellen, 
%       w = Vektor mit den zugehoerigen Gewichten
%function [x,w,count]=gaulob(n)
function [x,w]=gaulob(n)

m=floor((n-1)/2);
x=zeros(n,1);
w=zeros(n,1);
x(1)=1.0;
x(n)=-1.0;
w(1)=2./((n-1)*n);
w(n)=w(1);
%count=zeros(m,1);
for i=1:m
    z=cos(pi*i/(n-1));
    dz=1;
    while abs(dz)>eps    
        b=1;    
        c=z; 
        for k=2:n-1
            a=b;
            b=c;
            c=((2*k-1)*z*b-(k-1)*a)/k;
        end
        %dz=(z.*c-b)./(n*c+2*z*(z*c-b)/(1-z^2));
        dz=(z.*c-b)./(n*c);
        z=z-dz;   
        %count(i)=count(i)+1;         
    end
    x(i+1)=z;
    x(n-i)=-z;
    w(i+1)=2./((n-1)*n*c^2);
    w(n-i)=w(i+1);
end

end
