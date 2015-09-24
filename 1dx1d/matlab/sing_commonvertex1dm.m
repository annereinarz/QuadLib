%  Calculate Quadraturepoints for integrals of the form:
%
%  \int_{\tau_1}\int_{\tau_2} f(x,y) dxdy
%
%  where the function f has a singularity at [0,0]
%
%  INPUT:
%    n_gl  : the number of gauss legendre points
%    n_cgl : number for the composite gauss legendre quadrature
%               number of quadrature points is (n_cgl-1)*(n_cgl+2)/2+1
%
%  OUTPUT:
%  [X,W]   : Quadrature points on [0,1]^2 and quadrature weights
%
function [X,W]=sing_commonvertex1dm(n_gl, n_cgl)
if(n_gl==0 || n_cgl==0)
    error('need more than zero quadrature points');
end
[X1,W1]=gauleg(n_gl+1);
 X1=1/2+1/2*X1;
 W1=1/2*W1;
[X2,W2]=CGLquad(n_cgl); 
[X,W]=TensorQuad(X1',W1',X2,W2);
% Transformation of the quadrature points
W=W.*X(:,2);
X(:,1)=X(:,1).*X(:,2);

X = [X; X(:,2) X(:,1)];
W=[W; W];

return;
end