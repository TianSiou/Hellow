function [EW, EV, init_Vec, ritz_ew] = NullSp_SEP_Herm_JDSIRA_Driver_v2 (A, mtxdim, no_restart, ew_number, ...
                stop_tolerance, initial_vec, target, no_file, target_type, Orth_Number, CorEq, ...
                LSinfo, Precond_M_LS, Sigma_r, QBQ)

% mtx_A      : symmetric cofficient matrix of the eigenproblem                               
% mtxdim     : dimension of the cofficient matricx mtx_A      
% no_restart : number to restart the JD iteration        
% ew_number  : total number of eigenpairs wanted        

global zflag_ew flag_no_ew dim_eigsp

    % --- copy user define parameters

    RPPD                      = max(4, size(initial_vec,2)); 
    ProjProbDim               = size(initial_vec,2);
    init_Vec                  = zeros(mtxdim,ew_number+RPPD-1);
    init_Vec(:,1:ProjProbDim) = initial_vec; 
    
    flag_no_ew    = 0;
    dim_eigsp     = 0;
    EV            = zeros(mtxdim,ew_number); 

    RestartProjProbDim = RPPD;       % RestartProjProbDim must >= 1  

    zflag_ew           = zeros(ew_number+RPPD,1);

    % --- iteration for finding the first ew_number eigenpairs

    no_computed_ew = 1; 

    stop_flag = 1;
    
    while ( no_computed_ew <= ew_number && stop_flag ~= 0 ) 

          % --- compute the no_ite_ew-th eigenpair

          [no_computed_ew, EV, init_Vec, ProjProbDim, stop_flag, ritz_ew, n_root] = NullSp_SEP_Herm_JDSIRA_Locking_v2 ...
              ( A, EV, no_computed_ew, init_Vec, ProjProbDim, target, target_type, Orth_Number, ...
              no_file, RestartProjProbDim, no_restart, stop_tolerance, ew_number, CorEq, ...
              LSinfo, Precond_M_LS, Sigma_r, QBQ);
%           [no_computed_ew, EV, init_Vec, ProjProbDim, stop_flag, ritz_ew] = NullSp_SEP_Herm_JDSIRA_Locking_v2 ...
%               ( A, EV, no_computed_ew, init_Vec, ProjProbDim, target, target_type, ...
%               no_file, RestartProjProbDim, no_restart, stop_tolerance, ew_number, CorEq, ...
%               LSinfo, Precond_M_LS);
          
%           target = 0.9 * zflag_ew(no_computed_ew-1);

           %elapse_toc = etime(tarray)
           if ( ProjProbDim <= 0 ) 
               init_Vec(:,1) = randn(size(init_Vec,1),1); 
               rsdl          = norm(init_Vec(:,1));
               init_Vec(:,1) = init_Vec(:,1) / rsdl; 
               ProjProbDim   = 1;
           end
           
           %RestartProjProbDim     = max(4, RPPD-no_computed_ew+1);
           
%            if ( no_computed_ew > Orth_Number )
% %                target = 0.5 * max(real(zflag_ew(1:no_computed_ew-1,1)));
%                %target = 0.9 * real(zflag_ew(no_computed_ew-1,1)); %0.9
%                target = zflag_ew(no_computed_ew-5);
%            end
    if ( no_computed_ew > Orth_Number ) && mod(flag_no_ew  , Orth_Number) ==  0
%                        target = 0.5 * max(real(zflag_ew(1:no_computed_ew-1,1)));
        %target = 0.9 * real(zflag_ew(no_computed_ew-1,1)); %0.9
        target = zflag_ew(no_computed_ew-3);
    end
 
%     if n_root == 1
%        target = 0.9 * real(zflag_ew(no_computed_ew-1,1));
%     end

    end 
    
    init_Vec = init_Vec(:,1:ProjProbDim);
    
    EW = zflag_ew(1:no_computed_ew-1,1); %ew_number,1);
    fprintf(no_file,'kk = kk +1; \n\n');

