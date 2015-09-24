function [x,w]=CGLquad(n)

sigma = 0.1;    % parameter of geometric subdivision

x=[];
w=[];

xl = sigma;
xr = 1;

for j=1:n-1
    nj = n+(1-j);
    [x1,w1]=GLquad(nj,xl,xr);  
    
    x = [x; x1(1:nj)];
    w = [w; w1(1:nj)];
     
    xr=xl;
    xl=xl*sigma; 
end

nj = 1;
[x1,w1]=GLquad(nj,0,xr); 
x = [x; x1(1:nj)]; 
w = [w; w1(1:nj)];

if(abs(n*(n+1)/2-length(w))>eps), error('wrong length of quadrature points'); end

end
