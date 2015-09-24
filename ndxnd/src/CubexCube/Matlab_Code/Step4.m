function [t,wt]=Step4(N,k,t,wt)

for i=1:k
    if(N(i)~=0)
        t(:,N(i))=-t(:,N(i));
    end
end

end