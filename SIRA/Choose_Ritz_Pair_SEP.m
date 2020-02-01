function [ewno, ritz_ew, ritz_ev] = Choose_Ritz_Pair_SEP(iteno, ew, target, VRR, ...
                   Choose_Type, RestartProjProbDim, target_org, Choose_Type_org)

    % --- select the desired eigenpair (lambda,x) 
    %     such that lambda is the
    %     closest to target

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   New version, setup at Aug. 24, 2004
% 

warning off all

inVRR = VRR;
inew  = ew;
    ewp_idx     = zeros(iteno,1); 

    [sort_ritz, idx] = sort(ew);
%     VRR              = VRR(:,idx);

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
           dist_target(i) = abs(ew(i)-target);
           ewp_idx(i)     = i; 
       end
       jj = iteno+1;
    case('RGTR')
%      if the target value is 'Real' and the wanted eigenvalue is
%      'Greater' 'Than' target. The Ritz value in every iteration will 
%      be 'Real'
       for i = 1:iteno 
          if ( real(ew(i)) > real(target) && abs(imag(ew(i))) < 1.0e-10 ) 
              ewp_idx(jj)     = i;
              dist_target(jj) = abs(real(ew(i))-real(target));
              jj              = jj + 1; 
          end
       end
    case('RGTC')
%      if the target value is 'Real' and the real part of the wanted 
%      eigenvalue is 'Greater' 'Than' target. The Ritz value in every 
%      iteration can be 'Complex'
       for i = 1:iteno   
          if ( real(ew(i)) > real(target) )  
              ewp_idx(jj)     = i;
              dist_target(jj) = abs(ew(i)-target); 
              jj              = jj + 1; 
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
    
