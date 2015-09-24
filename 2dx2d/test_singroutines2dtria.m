function  test_singroutines2dtria

    make 2dtria; 
    
    epsilon=1/pi;
    
for counter=2
    switch(counter)
        case -1
            vertexlist=[1 0 1 -1 0 -1; 0 1 1 0 -1 -1];
            k=counter;  
        case 0
            vertexlist=[0 1 0 -1 0; 0 0 1 0 -1];
            k=counter; 
        case 1
            vertexlist=[0 0 1 -1; 0 1 0 0];
            k=counter;
        case 2
            vertexlist=[0 1 0; 0 0 1]; 
            k=counter;
    end
    
    % Test Convergence
    alpha = -4+k+epsilon;
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
    figure(1);
    plot(t(:,1),t(:,2),'r.'); hold on; plot(t(:,3),t(:,4), 'b.');
    sum(wt) 
    
    figure(3);
    semilogy(ndof.^(1/5), abs(1-Q/Q(end)),  '.-');
    hold all;
end

end
