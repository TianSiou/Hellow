function [ tau, beta, v] = householder_vector(alpha,x)
%  ZLARFG generates a complex elementary reflector H of order n, such
%  that
%
%        H' * ( alpha ) = ( beta ),   H' * H = I.
%             (   x   )   (   0  )
%
%  where alpha and beta are scalars, with beta real, and x is an
%  (n-1)-element complex vector. H is represented in the form
%
%        H = I - tau * ( 1 ) * ( 1 v' ) ,
%                      ( v )
%
%  where tau is a complex scalar and v is a complex (n-1)-element
%  vector. Note that H is not hermitian.
%
%  If the elements of x are all zero and alpha is real, then tau = 0
%  and H is taken to be the unit matrix.
%
%  Otherwise  1 <= real(tau) <= 2  and  abs(tau-1) <= 1 .
%

xnorm  = norm(x);
alphaR = real(alpha);
alphaI = imag(alpha);
if ( xnorm == 0 && alphaI == 0 ) 
    tau  = 0;
    v    = zeros(length(x),1);
    beta = alpha;
else
%
%   General case
%
    if ( alphaR == 0 )
        beta   = - dlapy3( alphaR, alphaI, xnorm );
    else
        beta   = - sign(alphaR) * dlapy3( alphaR, alphaI, xnorm );
    end
    tau    = ( beta - alphaR ) / beta - alphaI / beta * 1i;
    tmp    = alpha - beta;
    [p, q] = dladiv(1, 0, real(tmp), imag(tmp));
    alpha  = p + q * 1i;
    v      = alpha * x; 
end
% n     = length(x);
% v     = x;
% if ( alpha ~= 0 ) 
%     beta     = x(1) + sign(x(1)) * alpha;
%     v(2:n,1) = v(2:n,1) / beta;
% end
% v(1,1)  = 1; 
% beta    = 2 / ( 1 + (alpha^2 - x(1)^2)/beta^2 );
% alpha   = -sign(x(1)) * alpha;
