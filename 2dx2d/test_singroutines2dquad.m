function  test_singroutines2dquad

    make 2dquad; 

    epsilon=1/pi;
for counter=-1
    k=counter;
    switch(k)
        case -1
            vertexlist=[0 1 0  2 3 2;
                        0 0 1  0 0 1];            
        case 0
           vertexlist = [0  1 0 -1  0 ;
                         0  0 1  0 -1];
        case 1
            vertexlist=[0 1  0  0 ;
                        0 0  1 -1 ];
            
        case 2
            vertexlist=[0 0 1  ;
                        0 1 0.5];
    end
    
    % Test Convergence
    alpha = -2*2+k+epsilon;
    F = @(x,y) sqrt(sum((x-y).^2,2)).^alpha;
    Falt=@(z) sqrt(sum(z.^2,2)).^alpha;
    Falt2=@(x,y,z) sqrt(sum([z (x-y)].^2,2)).^alpha;

    cnt=1;
    for i=1:10
        if(k==-1)
            [t,wt]=squad2dquad(3*i,1, vertexlist);
        else
            [t,wt]=squad2dquad(i,2*i, vertexlist);
        end
        if(k==2)
            Q(cnt)=sum(Falt(t(:,5:6)).*wt);
        elseif(k==1)
            Q(cnt)=sum(Falt2(t(:,2), t(:,4), t(:,5)).*wt);
        else
            Q(cnt)=sum(F(t(:,1:2), t(:,3:4)).*wt);
        end
        ndof(cnt)=length(wt);
        cnt=cnt+1;
        
    end
    
     %Qer=convvalue(Q,ndof);
     %Qer=extrapolated_error(k);
 
      figure(2);
      semilogy(ndof.^(1/5), abs(1-Q/Q(end)),  '.-');
      hold all;
end
end

function a=convvalue(e,N)

un=e(end);
[b,c]=expconv(e(end-1:end), N(end-1: end), 5);
a=0.5*(un+c*exp(-b*N(end)^(1/5)));

end


function  [b,c]=expconv(e,N,gamma)

  % Find b
    b=-(log(e(2:end)./e(1:end-1)) ...
             ./(N(2:end).^(1/gamma)-N(1:end-1).^(1/gamma)));
  % Find c
    c= e(1:end-1)./exp(-b(1:end).*N(1:end-1).^(1/gamma));  
        
end

function a=extrapolated_error(k)

switch(k)
    case -1 
        a=0.014737672255389;
    case 0
        %a=0.805955657703852;
        a=1.811372822555820;
    case 1
        a=0.651945134328523;
    case 2
        a=8.994422639979659;
end

end