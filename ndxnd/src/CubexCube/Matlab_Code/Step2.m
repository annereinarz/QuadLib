function [zh,t,wt]=Step2(k,d,t,wt)

zh=t(:,1:k);
t(:,1:k)=zh+t(:,2*d-k+1 : 2*d);
t(:, d+1 : 2*d)=t(:,[(2*d-k+1 : 2*d)  (d+1 : 2*d-k)]);

end

