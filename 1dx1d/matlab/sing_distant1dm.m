%  Calculate Quadraturepoints for integrals of the form:
%
%  \int_{\tau_1}\int_{\tau_2} f(x,y) dxdy
%
%  where the function f has no singularity
%
%  INPUT:
%    n_gl  : the number of gauss legendre points
%    n_cgl : not applicable (unused in this type of quadrature)
%
%  OUTPUT:
%  [X,W]   : Quadrature points on [0,1]^2 and quadrature weights
%
function [X,W]=sing_none(n_gl, n_cgl)
if(n_gl==0 || n_cgl==0)
    error('need more than zero quadrature points');
end
[X1,W1]=gauleg(n_gl+1);
 X1=1/2+1/2*X1;
 W1=1/2*W1;
[X2,W2]=gauleg(n_gl+1);
 X2=1/2+1/2*X2;
 W2=1/2*W2;
[X,W] = TensorQuad(X1',W1',X2',W2');

end