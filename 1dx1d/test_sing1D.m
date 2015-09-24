function test_sing1D()

    make 1d;

    for counter=2:5
         switch(counter)
             case -1
                  a=0; b=1; c=-2; d=-1;
                  vertexlist=[a b  c d];
                  Qex=0.5*((b-d)^2*(3-log((b-d)^2)) ...
                            -(b-c)^2*(3-log((b-c)^2)) ...
                            -(a-d)^2*(3-log((a-d)^2)) ...
                            +(a-c)^2*(3-log((a-c)^2)));
                        k=-1;
            
              case 0
                   a=0; b=1; c=-1;
                   vertexlist=[a b c];
                   Qex=-0.5*((b-c)^2*(3-log((b-c)^2)) ...
                            -(a-c)^2*(3-log((a-c)^2)) ...
                            -(a-b)^2*(3-log((a-b)^2)));
                        k=0;

               case 1
                   a=0; b=1;
                   vertexlist=[a b];
                   Qex=-(a-b)^2*(3-log((a-b)^2));  
                   k=1;
             case 2
                 vertexlist=[0,1];
                 Qex=-35/72;
                 k=1;
             case 3
                 vertexlist=[0 1 -1];
                 Qex=1/72*(-29+48*log(2));
             case 4
                 vertexlist= [0 0; 1 0; 0 1]';
                 Qex=1/4*(-6+pi+log(4));
             case 5
                 vertexlist=[0 0; 1 0; sind(58) cosd(58)]';
                 Qex=-0.906076846593465;
         end
    
         % Test Convergence
         if(counter<2)
          F = @(x,y) log((x-y).^2);
         elseif(counter<4)
          F = @(x,y) log(abs(y-x)).*y.*y;
         else
             F = @(x,y) log(sqrt((y(:,1)-x(:,1)).^2+(y(:,2)-x(:,2)).^2));
         end
         
         cnt=1;
         for i=2:10
             [t,wt]=squad1d(0,i, 2*i,vertexlist);
             if(counter>=4)
                 Q(cnt)=sum(F(t(:,1:2), t(:,3:4)).*wt);
             else
                 Q(cnt)=sum(F(t(:,1), t(:,2)).*wt);
             end
             ndof(cnt)=length(wt);
             cnt=cnt+1;
         end
    
    
          %figure(1);
          %plot(t(:,1), t(:,2), '.');

          figure(2);
          semilogy(ndof.^(1/3), abs(1-Q/Qex),  '.-');
          hold all;
   end
end


