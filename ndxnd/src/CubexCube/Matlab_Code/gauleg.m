function [xg,wg] = gauleg(n)
%      [xg,wg] = gauleg(n)
%      find Gauss nodes xg(1),...xg(n) and weights wg(1),...,wg(n) for [-1,1]

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