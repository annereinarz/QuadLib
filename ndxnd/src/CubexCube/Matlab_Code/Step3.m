function [t,wt]=Step3(N,k,d,t,wt)

for i=1:k
    if(N(i)~=0)
        t(:,2*d-k+i)= t(:,N(i))+(1-abs(t(:,i))).* t(:,2*d-k+i);
    else
        t(:,2*d-k+i)=           (1-abs(t(:,i))).* t(:,2*d-k+i);
    end
    wt=(1-abs(t(:,i))).*wt;
end

end
