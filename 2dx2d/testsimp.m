function  testsimp

    make 2dtria; 
    
    epsilon=1/pi;
    Qh=0;
    
for counter=[221 211]
    switch(counter)
        case 221
            vertexlist=[0 1 1; 0 1 0];
            k=2;
        case 211
            vertexlist=[0 1 1 0; 0 1 0 1];
            k=1;
    end
    
    % Test Convergence
    alpha = -2+epsilon;
    F = @(z) sqrt(sum(z.^2,2)).^alpha; 
    
    cnt=1;
    for i=1:8
        [t,wt]=squad2dtria(i, 2*i, vertexlist);
        if(k==2)
            Q(cnt)=sum(F(t(:,5:6)).*wt);
        elseif(k==1)
            Q(cnt)=sum(F([t(:,5),t(:,4)-t(:,2)]).*wt);
        else
            Q(cnt)=sum(F(t(:,3:4) - t(:,1:2)).*wt);
        end
        ndof(cnt)=length(wt);
        cnt=cnt+1;
    end
    Qh=Qh+2*Q;
    Qh
end
    %figure(1);
    %plot(t(:,1),t(:,2),'r.'); hold on; plot(t(:,3),t(:,4), 'b.');
     
    figure(3);
    semilogy(ndof.^(1/5), abs(1-Qh/14.5558268793833356),  '.-');
    hold all;

end