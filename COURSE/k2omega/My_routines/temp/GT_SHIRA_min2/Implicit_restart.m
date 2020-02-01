function [new_H, new_R, new_Y, new_Z, unit_e] = Implicit_restart(converged, nmax, nstp, ...
                                                H, R, Y, Z)

% Get eigenvalues of Hessenberg matrix and sort by magnitude.
% Use smallest-magnitude eigenvalues as shifts. 

rtz                = eig(H(1:nmax,1:nmax), R(1:nmax,1:nmax));
[dummy,ndx]        = sort(abs(rtz));
rtz                = rtz(ndx);
new_Y              = Y;
new_Z              = Z;
nlast              = nmax - 2;
unit_e             = zeros(1,nmax);
unit_e(1, nmax)    = 1;
new_H              = H;
new_R              = R;
for jt = 1:nmax-nstp
   shift = rtz(jt); 
%
%  Do single shift QZ with shift rtz(jt)
%
   for j = converged:nlast   %converged = 0
      if j == converged
%   
%        Single shift M = H - a_1 * I
%        Compute first column of M 
%
         refvec = [ new_H(j+1,j+1) / new_R(j+1,j+1) - shift ;  new_H(j+2,j+1) / new_R(j+1,j+1) ];
      else
         refvec = new_H(j+1:j+2,j);
      end
%
%     Compute Givens rotation
%           | c       s |                  | x |     | r |
%       G = |           |   such that  G * |   |  =  |   |
%           |-conj(s) c |                  | y |     | 0 |
% 
%       where c is real, s is complex, and c^2 + |s|^2 = 1. 
%
%      G = givens(refvec(1), refvec(2)); 
      G = planerot(refvec);      
%
%     Update H:= G H  and  R:= G R
%
      if ( j == converged )
          new_H(j+1:j+2,j+1:nmax) = G * new_H(j+1:j+2,j+1:nmax);
          new_R(j+1:j+2,j+1:nmax) = G * new_R(j+1:j+2,j+1:nmax);
      else
          new_H(j+1:j+2,j:nmax) = G * new_H(j+1:j+2,j:nmax);
          new_H(j+2,    j     ) = 0;
          new_R(j+1:j+2,j:nmax) = G * new_R(j+1:j+2,j:nmax);
      end 
      
      refvec = new_R(j+2,j+1:j+2);
%      P      = givens(conj(refvec(1)), conj(refvec(2)));
      P      = planerot( [conj(refvec(1)); conj(refvec(2))] );
      P      = P';
      P      = [ P(:,2) P(:,1) ];
%
%     Update H := H * P, R := R * P   and    Y := Y * G', Z := Z * P
%   
      new_H (1:min(j+3,nmax),j+1:j+2) = new_H (1:min(j+3,nmax), j+1:j+2) * P; 
      new_R (1:j+2,          j+1:j+2) = new_R (1:j+2, j+1:j+2) * P;
      new_R (j+2,            j+1    ) = 0;
      new_Y (:,              j+1:j+2) = new_Y (:,     j+1:j+2) * G';
      new_Z (:,              j+1:j+2) = new_Z (:,     j+1:j+2) * P;
      unit_e(1,              j+1:j+2) = unit_e(1,     j+1:j+2) * P;
   end
   
end
