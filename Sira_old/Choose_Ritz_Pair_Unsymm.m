function [ewno, ritz_ew, ritz_ev] = Choose_Ritz_Pair_Unsymm(iteno, ew, ew_order, target, VRR, ...
                   vec_r_tol, Choose_Type, RestartProjProbDim, target_org, Choose_Type_org)

    % --- select the desired eigenpair (lambda,x) 
    %     such that lambda is the
    %     closest to target

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   New version, setup at Aug. 24, 2004
% 
global zflag_ew flag_no_ew

warning off all

    sort_ew     = zeros(iteno,1);
    ewp_idx     = zeros(iteno,1); 
    sort_No     = 5;

    sort_ritz = ew; 

    if ( ew_order > 1 )  
 
       mult_cn_ew    = zeros(flag_no_ew,1);  
       mult_ritz_w   = zeros(iteno,1);
       ld            = zeros(iteno,1); 
       idx_mult_ritz = zeros(iteno,1);

       sort_ew(1:flag_no_ew) = zflag_ew(1:flag_no_ew);

       % 
       % mult_cn_ew(i) : multiplicity of eigenvalue flag_ew(i)
       %
       % Find the multiplicity of eigenvalue. Suppose that
       %        \lambda_1 = ... = \lambda_k 
       % are mutiple eigenvalues. Then 
       %        mult_cn_ew(1) = k
       % and
       %        mult_cn_ew(j) = 0, for j = 2, ..., k
       %
       for i = 1:flag_no_ew
          if ( real(sort_ew(i)) > - 1.0e13 )  
             no_cn_ew  = 1;
             for i1 = i + 1:flag_no_ew
                if ( real(sort_ew(i1)) > - 1.0e13 )  
                   k        = -1;
                   flag_tol =  1;

                   while ( k <= 1 && flag_tol == 1 )
                      if ( abs(sort_ew(i1)-sort_ew(i)) <= 10.0^(k)*vec_r_tol  )  
                         flag_tol = 0;
                      else
                         k = k + 1;
                      end 
                   end

                   if ( flag_tol == 0 ) 
                      no_cn_ew    = no_cn_ew + 1;
                      sort_ew(i1) = -1.0e14;
                   end 
                end 
             end 
          else
             no_cn_ew = 0;
          end 
          mult_cn_ew(i) = no_cn_ew;
       end 

       % 
       % mult_ritz_w(i) : multiplicity of ritz value ALPHAR(i)
       % idx_ritz_w(j)  : denote the ritz value ALPHAR(j) to
       %                  be equal to ritz value ALPHAR(idx_ritz_w(j))
       %
       % Find the multiplicity of ritz value. Suppose that
       %        \theta_1 = ... = \theta_k 
       % are mutiple ritz values. Then 
       %        mult_ritz_w(1) = k
       % and
       %        mult_ritz_w(j) = 0, for j = 2, ..., k
       %
                            
       idx_ritz_w = zeros(iteno,1); 
       for i = 1:iteno
          if ( idx_ritz_w(i) == 0 ) 
             no_ritz_w = 1;
             for i1 = i + 1:iteno
                if ( idx_ritz_w(i1) == 0 )
                   k        = -1;
                   flag_tol =  1;
                   while ( k <= 1 && flag_tol == 1 )
                      if ( abs(sort_ritz(i1)-sort_ritz(i)) <= 10.0^(k)*vec_r_tol  )
                         flag_tol = 0;
                      else
                         k = k + 1;
                      end 
                   end
                   if ( flag_tol == 0 ) 
                      no_ritz_w      = no_ritz_w + 1;
                      idx_ritz_w(i1) = i;
                   end
                end
             end
          else
             no_ritz_w = 0;
          end
          mult_ritz_w(i) = no_ritz_w;
       end
       %
       % The ritz values which are equal to the convergent eigenvalues 
       % are set to be infty. 
       %
       for i = 1:flag_no_ew
          if ( real(sort_ew(i)) > - 1.0e13 ) 
             i1   = 1;
             flag = 1;
             while ( i1 <= iteno && flag == 1 )
                if ( idx_ritz_w(i1) == 0 ) 
                   k        = -1;
                   flag_tol = 1;    
                   while ( k <= sort_No && flag_tol == 1 )
                      if ( abs(sort_ritz(i1)-sort_ew(i)) <= 10.0^(k) * vec_r_tol  )
                         flag_tol = 0;
                         ind_sort = i1;
                      else
                         k = k + 1;
                      end
                   end
                   if ( flag_tol == 0 )
                      flag = 0;
                   else
                      i1 = i1 + 1;
                   end
                else
                   i1 = i1 + 1; 
                end
             end

             if ( flag == 1 )
                fprintf('Error in Choose_Ritz_Pair\n');
                fprintf('cn_ew(%3.0f) = \n', i, zflag_ew(i));
                fprintf('Convergent eigenvalues are:\n');
                disp(zflag_ew) 
                fprintf('\n');
                fprintf('Ritz values are:\n');
                disp(sort_ritz) 
                fprintf('\n');
             else
                if ( mult_cn_ew(i) == mult_ritz_w(ind_sort) ) 
                   sort_ritz(ind_sort) = 1.0e14;
                   for jj = 1:iteno
                      if ( idx_ritz_w(jj) == ind_sort )  
                         sort_ritz(jj) = 1.0e14;
                      end 
                   end
                   idx_ritz_w(ind_sort) = -1;
                elseif ( mult_cn_ew(i) < mult_ritz_w(ind_sort) ) 
                   %
                   % Another multiple eigenvalue which does not contain
                   % in the set of convergent eigenvalues belongs to
                   % the set of ritz values.
                   % The ritz values which are equal to the convergent 
                   % multiple eigenvalues are set to be infty so that
                   % the corresponding eigenvector does not be chosen. 
                   % How to choose these ritz values? Let span(V1) be the 
                   % eigenspace associated with the convergent eigenvalues
                   % and V = [ V1 V2 ], s = [ s1^T  s2^T ]^T. Then
                   %               u = V s
                   % is an approximated eigenvector. If s2 = 0, then u 
                   % is belong to span(V1) which is an already convergent
                   % eigenvector. Hence u is not a target vector. We want
                   % to choose u such that u is not belong to span(V1) and
                   % {V1, u} is linearly independent.
                   %
                   k                = 1;
                   idx_mult_ritz(k) = ind_sort;
                   ld(k)            = norm(VRR(ew_order:iteno,ind_sort));
                   for jj = 1:iteno
                      if ( idx_ritz_w(jj) == ind_sort ) 
                         k                = k + 1;
                         idx_mult_ritz(k) = jj;
                         ld(k)            = norm( VRR(ew_order:iteno,jj) );
                      end
                   end
                   if ( abs(k - mult_ritz_w(ind_sort)) > 0 ) 
                      fprintf('error in multiplicity of ritz value \n');
                   end
                   for jj = 1:mult_cn_ew(i)
                      tmp = ld(1);
                      k   = 1;
                      for i1 = 2:mult_ritz_w(ind_sort)
                         if ( tmp > ld(i1) )
                            tmp = ld(i1);
                            k   = i1;
                         end
                      end
                      sort_ritz(idx_mult_ritz(k)) = 1.0e14;
                      ld(k)                       = 1.0e14;
                   end
                   idx_ritz_w(ind_sort) = -1;
                else
                   sort_ritz(ind_sort) = 1.0e14;
                   for jj = 1:iteno
                      if ( idx_ritz_w(jj) == ind_sort )
                         sort_ritz(jj) = 1.0e14;
                      end
                   end
                   idx_ritz_w(ind_sort) = -1;
                   fprintf('Error for multiplicity of convergent \n');
                   fprintf('eigenvalue > multiplicity of ritz value \n');
                end
             end
          end
       end
    end

    %
    % Choose the ritz values which are closest the target value or 
    % approximate to the convergent eigenvalues.
    %
    jj          = 1;
    dist_target = zeros(iteno,1);

    switch Choose_Type
    case('CL')
%      choose the ritz value which is the closest target value
%
       for i = 1:iteno
          if ( abs(sort_ritz(i)) < 1.0e10 ) 
             ind_sort = 0;
             for k = 1:flag_no_ew
                if ( abs(sort_ritz(i)-zflag_ew(k)) < 1.0e3*vec_r_tol )  
                   ind_sort = k;
                end
             end
             if ( ind_sort > 0 )
                dist_target(jj) = abs(sort_ritz(i)-zflag_ew(ind_sort));
             else
                dist_target(jj) = abs(sort_ritz(i)-target);
             end
             ewp_idx(jj) = i;
             jj          = jj + 1;
          end
       end
    case('RGTR')
%      if the target value is 'Real' and the wanted eigenvalue is
%      'Greater' 'Than' target. The Ritz value in every iteration will 
%      be 'Real'
       for i = 1:iteno
          if ( abs(sort_ritz(i)) < 1.0e10 ) 
             if ( real(ew(i)) > real(target) && abs(imag(ew(i))) < 1.0e-12 ) 
                ewp_idx(jj)     = i;
                dist_target(jj) = abs(real(ew(i))-real(target));
                jj              = jj + 1;
             end
          end
       end
    case('RGTC')
%      if the target value is 'Real' and the real part of the wanted 
%      eigenvalue is 'Greater' 'Than' target. The Ritz value in every 
%      iteration can be 'Complex'
       for i = 1:iteno
          if ( abs(sort_ritz(i)) < 1.0e10 )  
             if ( real(ew(i)) > real(target) )  
                ewp_idx(jj)     = i;
                dist_target(jj) = abs(ew(i)-target); 
                jj              = jj + 1;
             end
          end
       end 
    otherwise 
       fprintf('In-correct Choose_Type; retype again \n');
    end 
  
    if ( jj == 1 ) 
       for i = 1:iteno
          ewp_idx(i)     = i;
          dist_target(i) = abs(sort_ritz(i)-target); 
       end
       jj = iteno + 1;
    end

    ewno = min(jj-1, RestartProjProbDim);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   New version, setup at Jan. 09, 2004
% 
    ritz_ev = zeros(iteno,ewno);
    ritz_ew = zeros(ewno,1);
    idx     = zeros(RestartProjProbDim,1);
    
    for i = 1:ewno
       ind_sort = 1;
       for k = 2:jj-1
          if ( dist_target(k) < dist_target(ind_sort) ) 
              if ( abs(target_org - target) < 1.0e-10 )
                  ind_sort = k;
              else
                  switch Choose_Type_org
                      case 'CL'
                          ind_sort = k;
                      case 'RGTR'
                          if ( real(ew(ewp_idx(k))) > real(target_org) && ...
                                  abs(imag(ew(ewp_idx(k)))) < 1.0e-12 )
                              ind_sort = k;
                          end
                      case 'RGTC'
                          if ( real(ew(ewp_idx(k))) > real(target_org) )
                              ind_sort = k;
                          end
                  end
              end
          end
       end

       if ( dist_target(ind_sort) > 1.0e12 )
           [~,ind_sort] = min( dist_target(1:jj-1) );
       end
       
       idx(i,1)              = ewp_idx(ind_sort);
       ritz_ew(i)            = ew(ewp_idx(ind_sort));
       ritz_ev(1:iteno,i)    = VRR(1:iteno,ewp_idx(ind_sort));
       ritz_ev(1:iteno,i)    = ritz_ev(1:iteno,i) / norm(ritz_ev(1:iteno,i));
       dist_target(ind_sort) = 1.0e14;
    end
 
%
%   modify at 2013/08/21
%
    if ( ewno < RestartProjProbDim )
        out_idx = zeros(iteno-ewno,1);
        kk      = 1;
        for ii = 1:iteno
            if ( isempty(find( abs(ii-idx) == 0, 1 )) )
                out_idx(kk) = ii;
                kk          = kk + 1;
            end    
        end
        
%      if the target value is 'Real' and the wanted eigenvalue is
%      'Greater' 'Than' target. The Ritz value in every iteration will 
%      be 'Real'
       jj = 1;
       for i = 1:iteno
          if ( abs(sort_ritz(i)) < 1.0e10 ) 
             if ( ~isempty(find( abs(i-out_idx) == 0, 1 )) )
                %if ( real(ew(i)) > real(target) ) 
                   ewp_idx(jj)     = i;
                   dist_target(jj) = abs(ew(i) - target);
                   jj              = jj + 1;
                %end
             end
          end
       end 
       
       [~,idx2] = sort(dist_target(1:jj-1));
       if ( ~isempty(idx2) )
          for ii = 1:min(RestartProjProbDim-ewno, length(idx2))
              ritz_ew(ewno+ii)         = ew(ewp_idx(idx2(ii)));
              ritz_ev(1:iteno,ewno+ii) = VRR(1:iteno,ewp_idx(idx2(ii)));
              ritz_ev(1:iteno,ewno+ii) = ritz_ev(1:iteno,ewno+ii) / ...
                  norm(ritz_ev(1:iteno,ewno+ii));
          end
          ewno = ewno + min(RestartProjProbDim-ewno, length(idx2));
       end
    end
    
