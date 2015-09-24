function main

% Input
type  = 21;
func_id  = 2;
epsilon  = 1/pi;

for i_type=1:length(type)
    cnt=1; 
    for nr=8
        switch type(i_type)
            case 1
                vertexlist=[0 1 0;
                            0 0 1]; 
                k=2;d=2;
            case 2
                vertexlist = [0  0   0.5;
                              0  0.5 0];
                k=2;d=2;
            case 3
                vertexlist=[ 0.5 0.5 0 1;
                             0   0.5 0 0];
                k=1;d=2;  
            case 4
                vertexlist=[0.5 1   0.5 0   0.5;
                            0.5 0.5 1   0.5 0  ];
                k=0;d=2;
           case 5
                 vertexlist=[ 0.5 0.7 0 1;
                         0   0.5 0 0.2];
                 k=1; d=2;
            case 19
                vertexlist = [0 1 0  -0.5 -1.5 -0.5  ;
                              0 0 1  -0.5 -0.5 -1.5 ];                
                d=2;  k=-1;
            case 20
                vertexlist = [0  1 0  -1  0  ;
                              0  0 1   0 -1 ];
                %vertexlist = [1 0 -1 4   5 ;
                %              1 2  0 0.5 2.5];
                d=2;  k=0;
            case 21
                vertexlist = [0 0  1  -1 
                              0 1  0   0 ];
                d=2;   k=1;
            case 22
                vertexlist = [0 0 1;
                              0 1 0];
                d=2;   k=2;
            case 29
                vertexlist = [0 0 1 0 -1 -2 -1 -1;
                              0 0 0 1 -1 -1 -2 -1;
                              0 1 0 0 -1 -1 -1 -2];
                d=3;   k=-1;
            case 30
                vertexlist = [0 0 1 0   0 -1  0;
                              0 0 0 1   0  0 -1;
                              0 1 0 0  -1  0  0];
                d=3;   k=0;
            case 31
                vertexlist = [0 1 0 0   0  0 ;
                              0 0 1 0  -1  0 ;
                              0 0 0 1   0 -1 ];
                d=3;   k=1;
            case 32
                vertexlist = [0 1 0 0  0 ;
                              0 0 1 0  0 ;
                              0 0 0 1 -1 ];
                d=3;   k=2;
            case 33
                vertexlist = [0 1 0 0;
                              0 0 1 0;
                              0 0 0 1];
                d=3;   k=3;
            case 39
                vertexlist = [0 1 0 0 0 -1 -2 -1 -1 -1;
                              0 0 1 0 0 -1 -1 -2 -1 -1;
                              0 0 0 1 0 -1 -1 -1 -2 -1;
                              0 0 0 0 1 -1 -1 -1 -1 -2];
                d=4;   k=-1;
            case 40
                vertexlist = [0 1 0 0 0  -1  0  0  0;
                              0 0 1 0 0   0 -1  0  0;
                              0 0 0 1 0   0  0 -1  0;
                              0 0 0 0 1   0  0  0 -1];
                d=4;   k=0;
            case 41
                vertexlist = [0 1 0 0 0  0  0  0;
                              0 0 1 0 0 -1  0  0;
                              0 0 0 1 0  0 -1  0 ;
                              0 0 0 0 1  0  0 -1 ];
                d=4;   k=1;
            case 42
                vertexlist = [0 1 0 0 0  0  0;
                              0 0 1 0 0  0  0;
                              0 0 0 1 0 -1  0;
                              0 0 0 0 1  0 -1];
                d=4;   k=2;
            case 43
                vertexlist = [0 1 0 0 0  0;
                              0 0 1 0 0  0;
                              0 0 0 1 0  0;
                              0 0 0 0 1 -1];
                d=4;   k=3;
           case 44
                vertexlist = [0 1 0 0 0;
                              0 0 1 0 0;
                              0 0 0 1 0;
                              0 0 0 0 1];
                d=4;   k=4;
            otherwise
                error('Unknown parallelotop location');
        end
        
        switch func_id
            case 1
               alpha = -2*d+k+epsilon;
               F = @(x,y,z,alpha) sqrt(sum((x-y).^2,2)).^alpha;
            case 5
               alpha = -d+epsilon;
               F = @(x,y,z,alpha) sqrt(sum((x-y).^2,2)).^alpha;
            case 2
                alpha = -2*d+k+epsilon;
                F = @(x,y,z,alpha) sqrt(sum(z.^2,2)).^alpha;
            case 11
                alpha = -2*d+k+epsilon;
                F = @(x,y,z,alpha) sqrt(sum(z.^2,2)).^-1;
            case 12
                alpha = -2*d+k+epsilon;
                F = @(x,y,z,alpha) sqrt(sum((x-y).^2,2)).^-1;
            otherwise
                error('Unknown func_id');
        end
        [Q(cnt)  ndof(cnt)]= squad(F,alpha,vertexlist,nr,2*nr,d,k);
        %[Q2(cnt) ndof2(cnt)]= squad_meth3(F,alpha,vertexlist,nr,2*nr,d,k);  %Only for k=d
        Q
        %Q2
         cnt=cnt+1;
    end
     %Qer=extrapolatedvalue(k,d);
     %Qer
     %figure(2);
     %semilogy(ndof.^(1/(2*d+1)), abs(1-Q/Qer),  '.-');
     %hold all;
     %semilogy(ndof2.^(1/(2*d+1)), abs(1-Q2/Qer),  '.-');
end

end

    function a=extrapolatedvalue(k,d)
        if(d==2)
            if(k==-1)
                a=0.050031938870694;
            elseif(k==0)
                a=1.811376905451605;
            elseif(k==1)
                a=5.8265665456293122;
            elseif(k==2)
                a=14.555627591816773;
            end
        elseif(d==3)
            if(k==-1)
                a=3.342063233068202e-04;
            elseif(k==0)
                a=0.466569689673609;
            elseif(k==1)
                a=2.080621269031285;
            elseif(k==2)
                a=8.331651803743892;
            elseif(k==3)
                a=28.397642748864797;
            end
        elseif(d==4)
            if(k==-1)
                a=1.0e-05*0.832703023021653;
            end
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