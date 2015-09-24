function  testsimp2

    make 2dtria; 
    
    epsilon=1/pi;
    Qh=0;
    Icoeff=[4 6 6];
    cnt=1;
    
for counter=[221 222 211 201]
    switch(counter)
        case 221
            vertexlist=[-1 0; 1 0; 0 sqrt(3)]';
            k=2;
        case 222
            vertexlist=[0 0; 1 0; 1/2 sqrt(3)/2]';
            k=2;
        case 211
            vertexlist=[0 0; 1/2 sqrt(3)/2; -1/2 sqrt(3)/2; 1 0]';
            k=1;
        case 201
            vertexlist=[1/2 sqrt(3)/2; 0 0; 1 0; -1/2 sqrt(3)/2; 0 sqrt(3)]';
            k=0;
    end
    
    % Test Convergence
    %alpha = -2+epsilon;
    %F = @(z) sqrt(sum(z.^2,2)).^alpha;
    F=@(z) log(sqrt(sum(z.^2,2)));
    
    if(counter==221)
        [t,wt]=squad2dtria(9, 18, vertexlist);
        Qtot=sum(F(t(:,5:6)).*wt);
    else       
        cnt2=1;
        for i=1:8
            [t,wt]=squad2dtria(i, 2*i, vertexlist);
            if(k==1)
                Q(cnt2)=sum(F([t(:,5),t(:,4)-t(:,2)]).*wt);
            else
                Q(cnt2)=sum(F(t(:,3:4) - t(:,1:2)).*wt);
            end
            ndof(cnt2)=length(wt);
            cnt2=cnt2+1;
        end
        Qh=Qh+Icoeff(cnt)*Q;
        cnt=cnt+1;
    end
end
    %figure(1);
    %plot(t(:,1),t(:,2),'r.'); hold on; plot(t(:,3),t(:,4), 'b.');
     
    %figure(3);
    semilogy(ndof.^(1/5), abs(1-Qh/Qtot(end)),  '.-');
    hold all;

end