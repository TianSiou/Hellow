function [x, no_iter, rsdl, flag] = newton_method_funhandle(name_fun_f, name_fun_df, x0, TOL)

narginchk(3,4);

if ( nargin < 4 )
    TOL = eps;
end

M = 1000 ;
k = 1;
x = x0 - name_fun_f(x0) / name_fun_df(x0) ;
while k <= M && abs(x-x0)/abs(x) >= TOL ;
    k  = k+1 ;
    x0 = x ;
    x  = x0 - name_fun_f(x0) / name_fun_df(x0) ;
end

if ( nargout >= 2 )
     no_iter = k;
     if ( nargout >= 3 )
         rsdl = name_fun_f(x);
         if ( nargout == 4 ) 
             if ( k <= M )
                 flag = 0; 
             else
                 flag = 1;
             end
         end
     end
end
  
% varargout(1) = {k};
% varargout(2) = {rsdl};
% if ( k <= M )
%     varargout(3) = {0};
% else
%     varargout(3) = {1};
% end