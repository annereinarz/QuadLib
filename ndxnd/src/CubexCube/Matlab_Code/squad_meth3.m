% this routine realizes a portion-wise computation and evaluation of the quadrature rule
function [Q ndof] = squad_meth3(F,alpha,vertexlist,nr,ns,d,k)

% parameters of the affine transformation to the physical coordinates
idx2 = [1:k+1 d+1+1:2*d-k+1];    % indices corresponding to the second cube
v10 = vertexlist(:,     1 );
v20 = vertexlist(:,idx2(1));
A1(:,1:d) = bsxfun(@minus,vertexlist(:,     2:d+1 ),v10);
A2(:,1:d) = bsxfun(@minus,vertexlist(:,idx2(2:d+1)),v20);

% compute quadrature rule on [0,1]x[-1,1]^d-1x[0,1]^d
[r,wr]=gauleg(nr);
[ralt,wralt]=gauleg(nr);
r=r'; wr=wr';ralt=ralt'; wralt=wralt';
for j=2:d-1
    [r,wr]=TensorQuad(ralt,wralt,r,wr);
end
[ralt,wralt]=GLquad(nr);
for j=d:2*d-1
    [r,wr] = TensorQuad(r,wr,ralt,wralt);
end
if(k==-1)
    [g,wg]=GLquad(nr); % only regular
else
    [g,wg]=CGLquad(ns);% one singular
end
[t,wt]=TensorQuad(g,wg,r,wr);

Q= 0; ndof=0;
% We assume k==d
  % Step 6: not necessary since check variables do not exist
  % Step 5: (r_1,tildep) in [0,1]x[-1,1]^k-1 to tildez in P_1
  [t(:,1:d),wt]=Step5_meth3(k, t(:,1:d),wt);
for j=1:d
    for sigma=1:2
      % Step 4: Permutations and Reflections
        [s, ws]=Step4_meth3(j,sigma,k,t,wt);
      % Step 3: tildeu in [0,1]^k to hatu in F_k(hat z)
        [s,ws]=Step3_meth3(k,d,s,ws);
      % Step 2: Transform J^kxJ^k to [-1,1]^k x F_k(zhat)
        [zh, s,ws]=Step2(k,d,s,ws);
      % Step 1: Transform to physical parallelotop
        [u,v,z,w] = Quad2PhyP(k,d,s(:,d+1:2*d),s(:,1:d)...
                        ,zh,ws,A1,v10,A2,v20);
        ndof=ndof+length(w);
        Q = Q + dot(F(u,v,z,alpha),w);
    end
end


end

