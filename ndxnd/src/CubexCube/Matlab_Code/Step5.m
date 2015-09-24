function [t,wt]=Step5(j,k,d,t,wt)

wt=(t(:,1).^(2*d-k-1)).*wt;
t(:,1:2*d-k)=bsxfun(@times, t(:,1), [t(:,2:j) ones(size(t(:,1))) t(:,j+1:2*d-k)]);

end