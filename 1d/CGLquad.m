function [x,w]=CGLquad(n)

sigma = 0.1;    % parameter of geometric subdivision
%b = 1.0;        % slope parameter
%delta = 1.0;    % Gevrey parameter

%m = ceil(b*n^(1/delta));

x=[];
w=[];

xl = sigma;
xr = 1;

for j=1:n-1
    %nj = ceil(n*(1+(1-j)/m)^delta);
    nj=n+1-j;     % Version without variable delta and b. The other version is not stable.
    [x1,w1]=gauleg(nj);
    x1=(xl+xr)/2+(xr-xl)/2*x1';
    w1=(xr-xl)/2*w1';
    
    x = [x; x1(1:nj)];
    w = [w; w1(1:nj)];
     
    xr=xl;
    xl=xl*sigma; 
end

%nj = ceil(n/(m^delta));
nj=n;
[x1,w1]=gauleg(nj);
 x1=xr/2+xr/2*x1';
 w1=xr/2*w1';   
x = [x; x1(1:nj)]; 
w = [w; w1(1:nj)];
end
