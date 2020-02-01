function sqrt_value = dlapy3(x, y, z)
%  DLAPY3 returns sqrt(x**2+y**2+z**2), taking care not to cause
%  unnecessary overflow.
%
      XABS = abs( x );
      YABS = abs( y );
      ZABS = abs( z );
      W    = max( [XABS YABS ZABS] );
      if ( W == 0 ) 
         sqrt_value = 0;
      else
         sqrt_value = W * sqrt( ( XABS / W )^2 + ( YABS / W )^2 + ( ZABS / W )^2 );
      end

