function main
make;
% Input
k=2; d=2;
func_id  = 2;

type  =  d*10+k;
for i_type=1:length(type)
    cnt=1; 
    for nr=1:5
         nr
         Q(cnt)= squad(nr,2*nr,type(i_type),func_id);   
         if(k==-1)
             ndof(cnt)=nr;
         else
             ndof(cnt)=(2^k*(2*d-k))^(1/(2*d+1))*nr;
         end
         
         cnt=cnt+1;
    end
     %Qer=convvalue(Q,ndof);
     Qer=extrapolatedvalue(k,d);
 
     figure(2);
     semilogy(ndof, abs(1-Q/Qer),  '.-');
     hold all;

end

end

    function a=extrapolatedvalue(k,d)
        if(d==2)
            if(k==-1)
                a=0.481333951582184; 
            elseif(k==0)
                a=0.748952218549046;
            elseif(k==1)
                a=1.112128689848246;
            elseif(k==2)
                a=2.973209598244779;
            end
        elseif(d==3)
            if(k==-1)
                a=0.288715087097229;
            elseif(k==0)
                a=0.578796997057372;
            elseif(k==1)
                a=0.708495100968968;
            elseif(k==2)
                a=0.980885087899876;
            elseif(k==3)
                a=1.882312098539536;
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