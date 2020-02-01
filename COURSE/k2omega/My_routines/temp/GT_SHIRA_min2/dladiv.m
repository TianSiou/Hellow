function [p, q] = dladiv(a, b, c, d)
%  DLADIV performs complex division in  real arithmetic
%
%                        a + i*b
%             p + i*q = ---------
%                        c + i*d
%

      if ( abs( d ) < abs( c ) ) 
         e = d / c;
         f = c + d * e;
         p = ( a + b * e ) / f;
         q = ( b - a * e ) / f;
      else
         e = c / d;
         f = d + c * e;
         p = ( b + a * e ) / f;
         q = ( - a + b * e ) / f;
      end

