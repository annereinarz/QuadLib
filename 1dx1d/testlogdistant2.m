function testlogdistant2

% Interval 1
x1=1;
% Interval 2
x2=1.3;
y2=0.3;
x3=0.8;
y3=1;

a=x2; b=-(x2-x3); c=y2; d=y3-y2; g=a-b*c/d;  h=b/d;

nr=17;
[t,wt]=GLquad(nr, 0, 1);
F=@(s)(atan((a+b*s)./(c+d*s)).*(y2+(y3-y2)*s));
Q= sum(F(t).*wt);

Ihh1=@(s)(0.5*s^2*atan(h+g/s)-g*(g-g*h^2)/(2*(h^2+1)^2)*atan(h+(h^2+1)/g*s)...
    -g^2*h/(2*(h^2+1)^2)*log((h^2+1)*s^2+2*g*h*s+g^2)+g/(2*(h^2+1))*s);
Ihh2=@(s)(s*atan(h+g/s)-2*g*h/(2*(h^2+1))*atan(h+(h^2+1)/g*s)+...
    (g/(2*(h^2+1)))*log(s^2+2*g*h/(h^2+1)*s+g^2/(h^2+1)));

Ih2=1/d^2*(Ihh1(y3)-Ihh1(y2))-y2/d^2*(Ihh2(y3)-Ihh2(y2));

Ih1=1/d*(Ihh2(y3)-Ihh2(y2));

I3=y2*Ih1+(y3-y2)*Ih2;

I3
Q

end
