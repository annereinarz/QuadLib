function [t,wt]=Step3_meth3(k,d,t,wt)

for i=1:k
    th=zeros(size(t(:,i)));
    I= t(:, i) <= 10^(-16);
    th(I)=t(I, i); 
    t(:,2*d-k+i)= -th+(1-abs(t(:,i))).* t(:,2*d-k+i);
    wt=(1-abs(t(:,i))).*wt;
end

end