function [z,w] = TensorQuad(x,u,y,v)
%computes tensor product of quadrature rules [x,u] and [y,v] 

[lx,nx] = size(x);	[ly,ny] = size(y);
    
z = zeros(lx*ly,nx+ny);
w = zeros(lx*ly,1);
% slow changing x
z(:,1:nx) = repmat(x,ly,1);
w         = repmat(u,ly,1);
% fast changing y
for i=1:ny
    y1        = repmat(y(:,i)',lx,1);
    z(:,nx+i) = reshape(y1,lx*ly,1);
end
v1 = repmat(v',lx,1);
w  = reshape(v1,lx*ly,1).*w;
end
