function [u,v,z,w] =Quad2PhyP(k,d,u,v,zh,w,A1,v10,A2,v20)

hat=1:k;
check=max(k,0)+1:d;

if(k<=0)
    z=(A2(:,check)*v(:,check)'-A1(:,check)*u(:,check)')';
elseif(k<d)
    z=zh*A1(:,hat)'+(A2(:,check)*v(:,check)')'-(A1(:,check)*u(:,check)')';
elseif(k==d)
    z=zh*A1(:,hat)';
end
if(k==-1)
  z=bsxfun(@plus, z, (v20-v10)');
end
u=bsxfun(@plus, v10, A1*u')';
v=bsxfun(@plus, v20, A2*v')'; 

w = w*abs(det(A1)*det(A2));

end