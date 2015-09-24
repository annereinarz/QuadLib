% this routine realizes a portion-wise computation and evaluation of the quadrature rule
function [Q ndof] = squad(F,alpha,vertexlist,nr,ns,d,k)

% parameters of the affine transformation to the physical coordinates
idx2 = [1:k+1 d+1+1:2*d-k+1];    % indices corresponding to the second cube
v10 = vertexlist(:,     1 );
v20 = vertexlist(:,idx2(1));
A1(:,1:d) = bsxfun(@minus,vertexlist(:,     2:d+1 ),v10);
A2(:,1:d) = bsxfun(@minus,vertexlist(:,idx2(2:d+1)),v20);

% compute quadrature rule on [0,1]^(2*d)
[r,wr]= GLquad(nr);
[ralt,wralt]=GLquad(nr);
for j=2:2*d-1
    [r,wr] = TensorQuad(ralt,wralt,r,wr);
end
if(k==-1)
    [g,wg]=GLquad(nr); % only regular
else
    [g,wg]=CGLquad(ns);% one singular
end
Q= 0; ndof=0;

[t,wt]=TensorQuad(g,wg,r,wr);
    
if(k==-1)
        % Transform to physical parallelotop
        [u,v,z,w] = Quad2PhyP(k,d,t(:,1:d),t(:,d+1:2*d),[],wt,A1,v10,A2,v20);
        Q = Q+dot(F(u,v,z,alpha),w);
        ndof=ndof+length(w);
elseif(k==0)
        for j=1:2*d
            %Step 5: Transform [0,1]^2d to D^N_j
            [s,ws]=Step5(j,k,d,t,wt);
            % Transform to physical parallelotop
            [u,v,z,w] = Quad2PhyP(k,d,s(:,1:d),s(:,d+1:2*d),[],ws,A1,v10,A2,v20);
            %figure(j);
            %plot(u(:,1), u(:,2), 'r.'); hold on; plot(v(:,1), v(:,2), 'b.');
            Q = Q + dot(F(u,v,z,alpha),w);
            ndof=ndof+length(w);
        end
else
        for numN=0:2^k-1
            N=FindSubset(numN,k);
            for j=1:2*d-k
                [s,ws]=BasicSubQuad(j,k,d,t,wt,N);
                %Step 3: utilde in[0,1]^k to uhat in F_N,M\N(zhat)
                [s,ws]=Step3(N,k,d,s,ws);
                %Step 2: Transform [-1,1]^k x F_k(zhat) to J^kxJ^k
                [zh, s,ws]=Step2(k,d,s,ws);
                %Step 1: Transform to physical parallelotop
                [u,v,z,w] = Quad2PhyP(k,d,s(:,d+1:2*d),s(:,1:d)...
                    ,zh,ws,A1,v10,A2,v20);
                ndof=ndof+length(w);
                Q = Q + dot(F(u,v,z,alpha),w);
                
            end
        end
end

end

