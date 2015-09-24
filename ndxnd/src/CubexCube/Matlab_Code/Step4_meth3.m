function [t,wt]=Step4_meth3(j,sigma,k,t,wt)
%Reflections
if(sigma==1)
    t(:,1)=-t(:,1);
end
%Permutations
th=t(:,1);
t(:,1)=t(:,j);
t(:,j)=th;

end