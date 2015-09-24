% Only works for k=d!
function [t,wt]=Step5_meth3(k,t,wt)

wt=t(:,1).^(k-1).*wt;
t(:,2:k)=bsxfun(@times, t(:,1), t(:,2:k));

end