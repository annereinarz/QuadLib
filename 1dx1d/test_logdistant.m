%
% For calculation of the exact solution see Documentation/log_distant.pdf
%
function test_logdistant

% Input the intervals
% Interval 1
x1=1;
% Interval 2
x2=1.3;
y2=0.3;
x3=0.8;
y3=1;

% Calculate r
r=(y3-y2)^2+(x3-x2)^2;

% Calculate m and n
m=-(1/r)*((x1-x2)*(x2-x3)+y2*(y3-y2));
n=(1/r)*((x1-x2)*(y3-y2)-y2*(x2-x3));

% Calculate e and f
e=-(1/r)*(x2*(x3-x2)+y2*(y3-y2));
f=(1/r)*(y2*x3-x2*y3);

if(y3==0 && y2==0)
      Ih=@(s) (s^2*log(s^2)-(s-1)^2*log((s-1)^2)+(1-2*s));
      I=  (sqrt(r))/(2*(x3-x2)) * (0.5*(Ih(x3)-Ih(x2)) )-sqrt(r);
else    
    if(n~=0)
        I1 = (x1-x2)* ( (1-m)*log((1-m)^2+n^2)-2+2*n*atan((1-m)/n) +m*log(m^2+n^2)+2*n*atan(m/n) )...
            -(x3-x2)*(0.5* (n^2-m^2+1)*log(n^2+(1-m)^2)- 0.5*(n^2- m^2)*log(n^2+m^2) +2*m*n* ...
            (atan((1-m)/n)+ atan(m/n))-0.5*(2*m+1))+(x1-0.5*(x2+x3))*log(r);
    else
        I1 = (x1-x2)*(log(r)-2*(m-1)*log(1-m)+2*m*log(-m)-2)...
            -(x3-x2)*(0.5*log(r) +0.5*(2*m^2*log(-m)-2*(m^2-1)*log(1-m)-2*m-1));
    end
    if(f~=0)
        I2 = x2* ( (1-e)*log((1-e)^2+f^2)-2+2*f*atan((1-e)/f) +e*log(e^2+f^2)+2*f*atan(e/f) )...
            +(x3-x2)*(0.5* (f^2-e^2+1)*log(f^2+(1-e)^2)- 0.5*(f^2- e^2)*log(f^2+e^2) +2*e*f* ...
            (atan((1-e)/f)+ atan(e/f))-0.5*(2*e+1))+(0.5*(x2+x3))*log(r);
    else
        I2 = x2*(log(r)-2*(e-1)*log(1-e)+2*e*log(-e)-2)...
            +(x3-x2)*(0.5*log(r) +0.5*(2*e^2*log(-e)-2*(e^2-1)*log(1-e)-2*e-1));
    end
    
    if(y3~=y2) % d!=0
       if(y2==0)  %c=0
          g=(x1-x2)/y3; h=(x2-x3)/y3;
          F2=@(s)(0.5*s^2*atan(h+g/s)-g*(g-g*h^2)/(2*(h^2+1)^2)*atan(h+(h^2+1)/g*s)...
                -g^2*h/(2*(h^2+1)^2)*log((h^2+1)*s^2+2*g*h*s+g^2)+g/(2*(h^2+1))*s);
          I3= 2*y3*(F2(1)-F2(0));
          
          g=x2/y3; h=(x3-x2)/y3;
          F2=@(s)(0.5*s^2*atan(h+g/s)-g*(g-g*h^2)/(2*(h^2+1)^2)*atan(h+(h^2+1)/g*s)...
                -g^2*h/(2*(h^2+1)^2)*log((h^2+1)*s^2+2*g*h*s+g^2)+g/(2*(h^2+1))*s);
          I4= 2*y3*(F2(1)-F2(0));
       else
           a=x1-x2; b=(x2-x3); c=y2; d=y3-y2; g=a-b*c/d;  h=b/d;

           Ihh1=@(s)(0.5*s^2*atan(h+g/s)-g*(g-g*h^2)/(2*(h^2+1)^2)*atan(h+(h^2+1)/g*s)...
               -g^2*h/(2*(h^2+1)^2)*log((h^2+1)*s^2+2*g*h*s+g^2)+g/(2*(h^2+1))*s);
           Ihh2=@(s)(s*atan(h+g/s)-2*g*h/(2*(h^2+1))*atan(h+(h^2+1)/g*s)+...
               (g/(2*(h^2+1)))*log(s^2+2*g*h/(h^2+1)*s+g^2/(h^2+1)));
           
           Ih2=1/d^2*(Ihh1(y3)-Ihh1(y2))-y2/d^2*(Ihh2(y3)-Ihh2(y2));
           
           Ih1=1/d*(Ihh2(y3)-Ihh2(y2));
           
           I3=y2*Ih1+(y3-y2)*Ih2;
           
           a=x2; b=-(x2-x3); c=y2; d=y3-y2; g=a-b*c/d;  h=b/d;

           Ihh1=@(s)(0.5*s^2*atan(h+g/s)-g*(g-g*h^2)/(2*(h^2+1)^2)*atan(h+(h^2+1)/g*s)...
               -g^2*h/(2*(h^2+1)^2)*log((h^2+1)*s^2+2*g*h*s+g^2)+g/(2*(h^2+1))*s);
           Ihh2=@(s)(s*atan(h+g/s)-2*g*h/(2*(h^2+1))*atan(h+(h^2+1)/g*s)+...
               (g/(2*(h^2+1)))*log(s^2+2*g*h/(h^2+1)*s+g^2/(h^2+1)));
           
           Ih2=1/d^2*(Ihh1(y3)-Ihh1(y2))-y2/d^2*(Ihh2(y3)-Ihh2(y2));
           
           Ih1=1/d*(Ihh2(y3)-Ihh2(y2));
           
           I4=y2*Ih1+(y3-y2)*Ih2;
       end
    elseif(y3==y2 && x2~=0)   %d=0, x2!=0
        a=x1-x2; b=x2-x3; c=y2;
        Ih1=@(s) (2*(a+b*s)*atan(abs((a+b*s)/c))-c*log(a^2+2*a*b*s+b^2*s^2+c^2))/(2*b);
        Ih2=@(s) (a*c*log(a^2+2*a*b*s+b^2*s^2+c^2)+...
            (a^2-c^2)*atan(abs(c/(a+b*s)))+b^2*s^2*atan(abs((a+b*s)/c))-b*c*s)/(2*b^2);
        I3 = 2*y2*(Ih1(1)-Ih1(0))+2*(y3-y2)*(1/(b^2))*(Ih2(1)-Ih2(0));
        
        a=x2; b=x3-x2; c=y2;
        Ih1=@(s)(2*(a+b*s)*atan((a+b*s)/c)-c*log(a^2+2*a*b*s+b^2*s^2+c^2))/(2*b);
        Ih2=@(s)(a*c*log(a^2+2*a*b*s+b^2*s^2+c^2)+...
            (a^2-c^2)*atan(c/(a+b*s))+b^2*s^2*atan((a+b*s)/c)-b*c*s)/(2*b^2);
        I4 = 2*y2*(Ih1(1)-Ih1(0))+2*(y3-y2)*(1/(b^2))*(Ih2(1)-Ih2(0));
    end
    I=sqrt(r)/2*(I1+I2+2*(I3+I4-x1));
end


% Approximate solution using tensor product quadrature
cnt=1;
for nr=2:10
    [r1,wr1]=GLquad(nr, 0, x1);   % quadrature on the first intervall
    [r2,wr2]=GLquad(nr, 0, 1);    % quadrature on the second intervall
    [t,wt]=TensorQuad(r1,wr1,r2,wr2);
    t=[t(:,1), zeros(size(t(:,1))), x2+(x3-x2)*t(:,2), y2+(y3-y2)*t(:,2)];
    wt=sqrt(r)*wt;
    plot(t(:,1), t(:,2), 'ro');
    hold on;
    plot(t(:,3),t(:,4),'bo');
    F=@(x,y)(log(sqrt((y(:,1)-x(:,1)).^2+(y(:,2)-x(:,2)).^2)));
    Q(cnt)=sum(F(t(:,1:2),t(:,3:4)).*wt);
    ndof(cnt)=nr^2;
    cnt=cnt+1;
end

I
Q

figure(2)
semilogy(ndof.^(1/2), abs(1-Q/I), 'r.-');


end
