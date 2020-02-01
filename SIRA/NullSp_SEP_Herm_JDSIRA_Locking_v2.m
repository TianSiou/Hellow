
function [ew_order, EV, init_Vec, ProjProbDim, stop_flag, ritz_ew, n_root] = NullSp_SEP_Herm_JDSIRA_Locking_v2(A, EV, ew_order, ...
    init_Vec, ProjProbDim, target, Choose_Type, Orth_Number, no_file, RestartProjProbDim, no_restart, ...
    vec_r_tol, ew_number, CorEq_org, LSinfo, Precond_M_LS , Sigma_r, QBQ)
%
% Use Jacobi-Davidson (JD) or Shift-Invert Residual Arnoldi (SIRA) Method to 
% compute the target eigenpair of the generalized eigenvalue problem
%                A x = lambda x
% where A is Hermitian.
%
global theta ev_crnt 
global zflag_ew flag_no_ew dim_eigsp
global ep_cpu_time ep_count sira_iter inner_count
% inner_count = 0;
if isa(A,'function_handle')
    MtxVecMult_A = A;
%     MtxVecMult_A = @(x) Sigma_r .* ( A(Sigma_r .* x) );
else
    % Check matrix and right hand side vector inputs have appropriate sizes
    [m,n] = size(A);
    if (m ~= n)
        fprintf('error:NonSquareMatrix\n');return;
    end
    MtxVecMult_A = @(x)Amtrix(x,A);
   
end

    % -------------------------
    % initialize matrix V, W, M
    % -------------------------
    
    tic;
    precond_flag       = LSinfo.precond;
    stop_flag          = 0; 
    EPIteNo            = 1;  
    iter               = 15;
    pcg_iter_number    = 1000; %1000;20
    prev_rsdl          = 1.0e2; 
    no_ew_near_tol     = 0;
    mtxdim             = size(init_Vec,1);
    n_root = 1;
    %
    % if the residual of ritz pair is small enough (less than tol_ritz_ew), then
    % we fix the target value (target_ritz_ew) as this ritz value.  
    % In the next iteration, we choose 
    % the ritz value which is closest to target_ritz_ew.
    %  
    % Modified at 2013/10/10
    %
    target_ritz_ew     = target;
    change_target      = 1;
    prev_ritz_ew       = 1.0e14;
    tol_ritz_ew        = 5.0e-3;
    Choose_Type_org    = Choose_Type;

    mtxV                  = zeros( mtxdim,       no_restart+2 ); 
    mtxM                  = zeros( no_restart+2, no_restart+2 );
    mtxV(:,1:ProjProbDim) = init_Vec(:,1:ProjProbDim);
    n_zero                = 1; 
%     Orth_cout             = 0;
%     if mod(flag_no_ew  , Orth_Number) ==  0 && flag_no_ew >= 2*Orth_Number
%         Orth_cout = flag_no_ew / Orth_Number -1;
%     else
%         Orth_cout             = 0;
%     end
    while ( EPIteNo <= 6000 && stop_flag == 0 )%%EPIteNo <= 6000
       
       start_time_ep = tic;
       %
       % Compute the Rayleigh quotient
       %      V^{*} A V   and   V^{*} B V
       %
       % MtxVecMult( V(:,ii), 0 ) : Compute A * V(:,ii)
       % MtxVecMult( V(:,ii), 1 ) : Compute B * V(:,ii)
       %
       for ii = 1:ProjProbDim
           mtx_A_prod_V              = MtxVecMult_A( mtxV(:,ii) );  
           mtxM(ii:ProjProbDim,ii)   = mtxV(:,ii:ProjProbDim)' * mtx_A_prod_V;
           mtxM(ii,ii+1:ProjProbDim) = mtxM(ii+1:ProjProbDim,ii)'; 
       end 

       % -----------------------------
       % iterate to find the eigenpair
       % -----------------------------

       iteno      = ProjProbDim; 
       while (iteno <= no_restart && stop_flag == 0 ) 

          % --- solve for the expanded generalized eigensystem
          
          [VRR, ew] = eig(mtxM(1:iteno, 1:iteno));
          ew        = diag(ew);

          % --- select the desired eigenpair (lambda,x) 
          %     where lambda is the
          %     closest to (real_target, imag_target)

          [ewno, ritz_ew, ritz_ev] = Choose_Ritz_Pair_SEP(iteno, ew, ...
              target_ritz_ew, VRR, Choose_Type, RestartProjProbDim, ...
              target, Choose_Type_org);  
%           norm(mtxM(1:iteno, 1:iteno)*ritz_ev(:,1) - ritz_ew(1) *ritz_ev(:,1));
          idx_ew      = find(abs(ritz_ew - ritz_ew(1)) < 1.0e-6);
          n_zero_prev = n_zero;
          n_zero      = max(1,length(idx_ew));
          
          % ----------------------------------------------------
          % comput u (eigenvector of the PEP), p = A'(labda)u, 
          %        r = A(lambda)u (residual of the PEP)
          % ----------------------------------------------------

              ev_crnt_tmp = mtxV(:,1:iteno) * ritz_ev(1:iteno,idx_ew(1:n_zero)); 
              for ii = 1:n_zero
                  ev_crnt_tmp(:,ii) = ev_crnt_tmp(:,ii) / norm(ev_crnt_tmp(:,ii));
              end
%               rsdl_vec_tmp = MtxVecMult_A(ev_crnt_tmp);
%               rsdl_vec_tmp = rsdl_vec_tmp - ev_crnt_tmp * diag(ritz_ew(idx_ew(1:n_zero)));
              for p = 1 : n_zero        %%%modify by siu
                  rsdl_vec_tmp(:,p) = MtxVecMult_A(ev_crnt_tmp(:,p));
                  rsdl_vec_tmp(:,p) = rsdl_vec_tmp(:,p) - ev_crnt_tmp(:,p) * ritz_ew(p);
              end
              rsdl_tmp     = zeros(n_zero,1);
              for ii = 1:n_zero
                  rsdl_tmp(ii,1) = norm(rsdl_vec_tmp(:,ii));
              end
              
              if ( n_zero > 1 )
                  [~,idx_rsdl] = sort(rsdl_tmp);
                  rsdl_tmp     = rsdl_tmp(idx_rsdl,1);
                  rsdl_vec_tmp = rsdl_vec_tmp(:,idx_rsdl);
                  ev_crnt_tmp  = ev_crnt_tmp(:,idx_rsdl);
                  
                  ritz_ew(idx_ew(1:n_zero))         = ritz_ew(idx_ew(idx_rsdl));
                  ritz_ev(1:iteno,idx_ew(1:n_zero)) = ritz_ev(1:iteno,idx_ew(idx_rsdl));
                      
              end
%               Rsdl    = [Rsdl, rsdl_tmp];
              
              print_result = 1;
              for ii = 1:n_zero    
                  if ( print_result == 1 )
                      if ( imag(ritz_ew(idx_ew(ii))) > 0.0  ) 
                          fprintf('ritz_ew(%3.0f, %2.0f,%4.0f) = %18.12e+%18.12ei; rsdl(%2.0f,%4.0f) = %8.2e;\n', ...
                              ew_order, ii, EPIteNo, real(ritz_ew(idx_ew(ii))), imag(ritz_ew(idx_ew(ii))), ...
                              ew_order, EPIteNo, rsdl_tmp(ii,1));
                      elseif ( imag(ritz_ew(idx_ew(ii))) < 0  )
                          fprintf('ritz_ew(%3.0f, %2.0f,%4.0f) = %18.12e-%18.12ei; rsdl(%2.0f,%4.0f) = %8.2e;\n', ...
                              ew_order, ii,EPIteNo,real(ritz_ew(idx_ew(ii))),abs(imag(ritz_ew(idx_ew(ii)))), ...
                              ew_order, EPIteNo, rsdl_tmp(ii,1));
                      else
                          fprintf('ritz_ew(%3.0f, %2.0f,%4.0f) = %18.12e; rsdl(%2.0f,%4.0f) = %8.2e;\n', ...
                              ew_order, ii, EPIteNo, real(ritz_ew(idx_ew(ii))), ...
                              ew_order, EPIteNo, rsdl_tmp(ii,1));
                      end 
                  end
          
              end 
              
          ev_crnt  = ev_crnt_tmp(:,1);
          rsdl_vec = rsdl_vec_tmp(:,1);
          rsdl     = rsdl_tmp(1,1);
          
%           % =====
%           if ( n_zero > 1 )
%           tt = mtxV(:,1:iteno) * ritz_ev(1:iteno,1);
%           tt = tt / norm(tt);
%           fprintf('Check error of ritz vector = %11.4e \n',norm(tt-ev_crnt));
%           fprintf('Check new residual = %11.4e, rsdl_org = %11.4e \n', norm(MtxVecMult_A(tt) - ritz_ew(1) * tt), rsdl);
%           end
%           % =====
          
%           %         
%           %  ev_crnt : current target approximated eigenvector
%           %
%           ev_crnt  = mtxV(:,1:iteno) * ritz_ev(1:iteno,1);  
%           ev_crnt  = ev_crnt / norm( ev_crnt); 
%           %
% 	  % Compute the residual of the ritz pair ( ritz_ew(1), ev_crnt )
% 	  %
%           rsdl_vec = MtxVecMult_A(ev_crnt);
%           rsdl_vec = rsdl_vec - ritz_ew(1) * ev_crnt; 
%           rsdl     = norm(rsdl_vec);
%   
%           print_result = 1;
%           if ( print_result == 1 )
%              if ( imag(ritz_ew(1)) > 0.0  ) 
%                 fprintf('ew_crnt(%2.0f,%4.0f) = %18.12e+%18.12ei; rsdl(%2.0f,%4.0f) = %8.2e;\n', ...
%                    ew_order, EPIteNo, real(ritz_ew(1)), imag(ritz_ew(1)), ...
%                    ew_order, EPIteNo, rsdl);
%              elseif ( imag(ritz_ew(1)) < 0  )
%                 fprintf('ew_crnt(%2.0f,%4.0f) = %18.12e-%18.12ei; rsdl(%2.0f,%4.0f) = %8.2e;\n', ...
%                    ew_order,EPIteNo,real(ritz_ew(1)),abs(imag(ritz_ew(1))), ...
%                    ew_order, EPIteNo, rsdl);
%              else
%                 fprintf('ew_crnt(%2.0f,%4.0f) = %18.12e; rsdl(%2.0f,%4.0f) = %8.2e;\n', ...
%                    ew_order, EPIteNo, real(ritz_ew(1)), ...
%                    ew_order, EPIteNo, rsdl);
%              end 
%           end
          
          if ( rsdl <= 10.0e0 * vec_r_tol && rsdl > vec_r_tol )  
             no_ew_near_tol = no_ew_near_tol + 1;
          else
             no_ew_near_tol = 0;
          end

          % -----------------
          % check convergence
          % -----------------
          
          if ( rsdl <= vec_r_tol || no_ew_near_tol >= 5 ) 
             stop_flag = 1; 
             n_zero    = n_zero_prev;
          else  
             if ( rsdl < tol_ritz_ew && change_target == 1 && ...
                  abs(prev_ritz_ew - ritz_ew(1))/abs(ritz_ew(1)) < 5.0e-2 )

                target_ritz_ew = ritz_ew(1);
                change_target  = 0;
                Choose_Type    = 'CL';
             end
             prev_ritz_ew = ritz_ew(1);
%
%            modify at 2007, 1, 19
%
%             if ( iteno + increment <= no_restart ) then
%              if ( iteno + 1 <= dim_eigsp+no_restart ) 
             if ( iteno + n_zero <= no_restart )
%
%               Compute the appending vector
%
                %
                % Adaptive determine the iteration numbers for solving correction 
                % equation in JD method
                %

                if ( rsdl >= 1.0e-1 && EPIteNo < 14 ) 
                   iter            = 50;
%                    iter_flag       = 1;
                   pcg_iter_number = 50; %10^3;
%                 elseif ( iter_flag == 1 )  
%                    iter            = 30; %15;
%                    iter_flag       = 0;
%                    pcg_iter_number = 100; %1000;
                elseif ( rsdl < 1.0e-1 && prev_rsdl/rsdl < 4.0e0 ) 
                   iter            = min(500, iter+10); %min(1000, iter+20);
                   pcg_iter_number = min(500,pcg_iter_number+10);
                end 
                prev_rsdl = rsdl;

%
%                  init_Vec(:,1) := M_A^{-1} r = M_A^{-1} tmpVt(:,2)
%                  tmpVt(:,2)    := M_A^{-1} p = M_A^{-1} tmpVt(:,1)
%
                theta = ritz_ew(1);
                
                if ( ew_order == 1 && EPIteNo <= 10 && rsdl > 5.0e-3 )
                    CorEq = 'SIRA';%'JD'; %'SIRA';
                    LSinfo.precond = 'no';
                    if ( strcmp(CorEq_org,'JD') || strcmp(CorEq_org,'jd') )
                        LSinfo.precond = 'no';
                    end
                else
                    CorEq          = CorEq_org; %'JD'; %'SIRA';
                    LSinfo.precond = precond_flag;
                    %LS_solver      = LSinfo.solver;
                end
                %CorEq = CorEq_org;
                
                switch CorEq
                    case 'JD'
                        LS_solver     = LSinfo.solver;
                        init_Vec(:,1) = solve_approx_vec_t_SEP(A, Precond_M_LS, LS_solver, theta, ...
                            rsdl_vec, ev_crnt, ew_order, iter, rsdl); 
                           
                    case 'SIRA' 

                        linear_system_tol = 1.0e-3; %%1.0e-3
                        %pcg_iter_number   = 100; %iter; %2 * 10^3;
                        LS_solver         = LSinfo.solver; %'minres';
                                                  
                        if ( strcmp(LSinfo.precond, 'yes') || strcmp(LSinfo.precond, 'YES') )
                            init_Vec(:,1) = SIRA_solve_residual_LS_SEP(A, target, rsdl_vec, ...
                                LS_solver, pcg_iter_number, linear_system_tol, ev_crnt, Precond_M_LS); 
                        else
%                             init_Vec(:,1) = SIRA_solve_residual_LS_SEP(A, target, rsdl_vec, ...
%                                 LS_solver, pcg_iter_number, linear_system_tol); %, ev_crnt); 
                             init_Vec(:,1) = SIRA_solve_residual_LS_SEP(QBQ, target, rsdl_vec, ...
                                LS_solver, pcg_iter_number, linear_system_tol, Sigma_r);      
                        end
                          n_root = n_zero; %%modfy by siu
                          n_zero = 1;   
                          for ii = 2:n_zero
                            if ( strcmp(LSinfo.precond, 'yes') || strcmp(LSinfo.precond, 'YES') )
                                init_Vec(:,ii) = SIRA_solve_residual_LS_SEP(A, target, rsdl_vec_tmp(:,ii-1), ...
                                    LS_solver, pcg_iter_number, linear_system_tol, ev_crnt, Precond_M_LS); 
                            else
%                                 init_Vec(:,ii) = SIRA_solve_residual_LS_SEP(A, target, rsdl_vec_tmp(:,ii-1), ...
%                                     LS_solver, pcg_iter_number, linear_system_tol); %, ev_crnt); 
                                init_Vec(:,ii) = SIRA_solve_residual_LS_SEP(QBQ, target, rsdl_vec_tmp(:,ii), ...   %%modify by siu
                                    LS_solver, pcg_iter_number, linear_system_tol, Sigma_r);
                            end
                          end
                end

%                 if ( n_zero > 1 )
%                     [i, init_Vec(:,1:i)] = Gram_Schmidt(n_zero, init_Vec(:, 1:n_zero));
%                     %[mtx_Q,~] = qr(init_Vec(:, 1:n_zero),0); 
%                     %init_Vec(:, 1:n_zero) = mtx_Q;
% %                     if ( n_zero ~= i )
% %                         fprintf('not full rank of corrected vectors \n');
% %                     end
%                 end
                    % ----------------------------------------------
                    % orthogonalize t and then nomalize it
                    % ----------------------------------------------
                    if(n_zero > 1)
                       for k = 2 : n_zero
                           init_Vec(:,k) = init_Vec(:,k) - ((init_Vec(:,k)' * init_Vec(:, 1 : k-1)) * init_Vec(:, 1 : k-1)' )';
                       end
                    end
                for ii = 1:n_zero
                    rsdl          = norm(init_Vec(:,ii));
                    if ( rsdl < 1.0e-9 )  
                        init_Vec(:,ii) = init_Vec(:,ii) / rsdl;
                    end
%                 for ii = 1:1
%                     rsdl          = norm(init_Vec(:,ii));
%                     if ( rsdl < 1.0e-9 )  
%                         init_Vec(:,ii) = init_Vec(:,ii) / rsdl;
%                     end
                    
                    % ----------------------------------------------
                    % orthogonalize t against V and then nomalize it
                    % ----------------------------------------------
%                     if ( ew_order > 1 )
% %                         for jj = 1:ew_order-1
% %                             ew_tmp         = EV(:,jj)' * init_Vec(:,ii);
% %                             init_Vec(:,ii) = init_Vec(:,ii) - ew_tmp * EV(:,jj);
% %                         end
%                         init_Vec(:,ii) = init_Vec(:,ii) - ((init_Vec(:,ii)' * EV(:, 1 : ew_order-1)) * EV(:, 1 : ew_order-1)' )';
%                     end
                    if flag_no_ew > 1 && flag_no_ew < Orth_Number
                    init_Vec(:,ii) = init_Vec(:,ii) - ((init_Vec(:,ii)' * EV(:, 1 : flag_no_ew)) * EV(:, 1 : flag_no_ew)' )';
                    elseif flag_no_ew >= Orth_Number
                        init_Vec(:,ii) = init_Vec(:,ii) - ( (init_Vec(:,ii)' * EV(:, flag_no_ew - Orth_Number + 1 : flag_no_ew)) * EV(:, flag_no_ew - Orth_Number + 1 : flag_no_ew)' )';
                    end

%                     if flag_no_ew >= Orth_Number
%                             init_Vec(:,ii) = init_Vec(:,ii) - (init_Vec(:,ii)' * EV(:, Orth_cout*Orth_Number + 1 : (Orth_cout+1)*Orth_Number ) * EV(:, Orth_cout*Orth_Number+1 : (Orth_cout+1)*Orth_Number)' )';
% %                         if mod(flag_no_ew  , Orth_Number) ==  0 && flag_no_ew >= 2*Orth_Number
% %                            Orth_cout = Orth_cout + 1;
% %                         end
%                     end
%                     for j = 1:iteno+ii-1
%                         ew_tmp         = mtxV(:,j)' * init_Vec(:,ii);
%                         init_Vec(:,ii) = init_Vec(:,ii) - ew_tmp * mtxV(:,j);
%                     end
                    init_Vec(:,ii) = init_Vec(:,ii) - ((init_Vec(:,ii)' * mtxV(:, 1 : iteno+ii-1)) * mtxV(:, 1 : iteno+ii-1)' )';
                    
                    mtxV(:,iteno+ii) = init_Vec(:,ii) / norm(init_Vec(:,ii));
                    if norm( mtxV(:, iteno+ii) - mtxV(:, iteno)) < 1e-1
                        a = 0;
                    end
                    if abs(mtxV(:,iteno+ii)'  * mtxV(:,iteno) ) > 1e-6 %%check v is not converge last v
                        a = 0;
                    end
%                     abs(mtxV(:,iteno+ii)'  * mtxV(:,iteno) )
                    %%% iteno + ii  is k+1
                    % ---------------
                    % update V^T A V  and  V^T B V
                    % --------------- 
                    mtx_A_prod_V                = MtxVecMult_A( mtxV(:,iteno+ii) ); 
%                     mtx_A_prod_V                = MtxVecMult_A( mtxV(:,1:iteno+ii) ); 
                    mtxM(1:iteno+ii,iteno+ii)   = mtxV(:,1:iteno+ii)' * mtx_A_prod_V;
                    mtxM(iteno+ii,1:iteno+ii-1) = mtxM(1:iteno+ii-1,iteno+ii)';
                    
                end
             end 

          end 

%           iteno   = iteno + n_zero;
          iteno   = iteno + 1;
          EPIteNo = EPIteNo + 1;

       end

       
       % --- construct the initial search space with ProjProbDim Ritz pairs 
       %     for restarting
       %
       % The orthonormal basis for the convergent eigenvactors is always locked at 
       % the first "dim_eigsp" columns of the matrix mtxV. The current approximate 
       % eigenvector is put on the "dim_eigsp"+1-th column of mtxV. 
       % We also collect the Ritz vectors for the Ritz values which are the second, 
       % third, ... closest to the target value to accelerate the convergence in computing	
       % next target eigenpair.
       %
       if ( stop_flag == 0 )

          mtxM  = zeros( no_restart+2, no_restart+2 );
          iteno = iteno - n_zero;

          jj               = min(RestartProjProbDim-1,ew_number+RestartProjProbDim-1);
          init_Vec(:,1:jj) = mtxV(:,1:iteno) * ritz_ev(1:iteno,2:jj+1);

%           k             = dim_eigsp; 
%           mtxV(:, k+1)  = ev_crnt;
%           k             = k + 1; 
% 
%           mtxV(:,k+1:k+jj) = init_Vec(:,1:jj);
%           k                = k + jj;
% 
% %           i            = dim_eigsp; 
%           [i, mtxV(:,1:i)] = Gram_Schmidt(k, mtxV(:, 1:k));
%           ProjProbDim      = i;
          % ==
          if ( ew_order > 1 )
              mtxV(:,1:ew_order-1) = EV(:,1:ew_order-1);
          end
          mtxV(:, ew_order)              = ev_crnt;
          mtxV(:,ew_order+1:ew_order+jj) = init_Vec(:,1:jj);
          jj                             = jj + ew_order;

          [i, mtxV(:,1:i)]      = Gram_Schmidt(jj, mtxV(:, 1:jj));
          ProjProbDim           = i - ew_order + 1;
          mtxV(:,1:ProjProbDim) = mtxV(:,ew_order:i);
          

       end 
       ep_cpu_time(ep_count) = toc(start_time_ep);
    end 
           sira_iter(ep_count)   = EPIteNo;
           if length(sira_iter) ~= ep_count
               a = 0;
           end
    %
    %   Compute Ritz pairs and the corresponding residual
    %
    rsdl_check = 1;
    tmp_rsdl   = ones(ewno,1); 

    jj         = 0;
    if n_zero > 1
        a = 0;
    end  
    for i = 2:ewno %RestartProjProbDim
          if ( jj+1 <= ew_number+RestartProjProbDim-1 )   
             init_Vec(:, i) = mtxV(:,1:iteno-n_zero) * ritz_ev(1:iteno-n_zero,i);
%               init_Vec(:, i) = mtxV(:,1:iteno-1) * ritz_ev(1:iteno-1,i);  %% modify by siu
             init_Vec(:, i) = init_Vec(:,i) / norm(init_Vec(:,i)); 
             jj                = jj + 1;
             if ( rsdl_check == 1 )   
                 rsdl_vec   = MtxVecMult_A(init_Vec(:,i)) - ritz_ew(i) * init_Vec(:,i);
                tmp_rsdl(i) = norm(rsdl_vec); 
                
                if ( imag(ritz_ew(i)) > 0 ) 
                   fprintf('rsdl(%24.16e+%24.16ei) = %14.4e \n', real(ritz_ew(i)), ...
                       imag(ritz_ew(i)),tmp_rsdl(i));
                elseif ( imag(ritz_ew(i)) < 0 )
                   fprintf('rsdl(%24.16e-%24.16ei) = %14.4e \n', real(ritz_ew(i)), ...
                       -imag(ritz_ew(i)),tmp_rsdl(i));
                else
                   fprintf('rsdl(%24.16e) = %14.4e \n', real(ritz_ew(i)),tmp_rsdl(i));
                end
             end
          end
    end
    
    flag_no_ew = flag_no_ew + 1;
    k          = dim_eigsp;

    %
    % Save the converging eigenpair and lock it on matrix V
    %

       tmpVt(:,1)           = ev_crnt / norm(ev_crnt);
       tmpVt(:,2)           = zeros(size(ev_crnt,1),1);
       if ( rsdl_check == 1 ) 
           rsdl_vec = MtxVecMult_A(tmpVt(:,1)) - ritz_ew(1) * tmpVt(:,1);
           rsdl     = norm(rsdl_vec);
           
           if ( imag(ritz_ew(1)) > 0 ) 
               fprintf('rsdl(%24.16e+%24.16ei) = %14.4e \n', real(ritz_ew(1)), ...
                       imag(ritz_ew(1)),rsdl);
           elseif ( imag(ritz_ew(1)) < 0 )
               fprintf('rsdl(%24.16e-%24.16ei) = %14.4e \n', real(ritz_ew(1)), ...
                       -imag(ritz_ew(1)),rsdl);
           else
               fprintf('rsdl(%24.16e) = %14.4e \n', real(ritz_ew(1)),rsdl);
           end
                
       end 
%        mtxV(:,k+1  )        = tmpVt(:,1);
       EV(:,k+1)            = tmpVt(:,1);
       k                    = k + 1;
       zflag_ew(flag_no_ew) = ritz_ew(1);
       dim_eigsp            = dim_eigsp + 1;
    %
    %   Check any other converging eigenpair. If yes, lock it again into V
    %
    flag_cn = zeros(ewno,1);
    j1      = 0;
    j       = 0;
    for i = 2:ewno %RestartProjProbDim 
          if ( j1+1 <= ew_number+RestartProjProbDim-1 ) 
             if ( tmp_rsdl(i) <= vec_r_tol)     %%%both two ew converge
%                 mtxV(:,k+1)          = init_Vec(:,j1+1);
                EV(:,k+1)            = init_Vec(:,i);
                k                    = k + 1;
                flag_cn(i)           = 1;
                j                    = j + 1;
                flag_no_ew           = flag_no_ew + 1;
                zflag_ew(flag_no_ew) = ritz_ew(i);
                dim_eigsp            = dim_eigsp + 1;
                j1                   = j1 + 1;
                output(no_file, ew_order+j, 0, 0.0d0, ritz_ew(i), tmp_rsdl(i));
             end
          end 
    end

    output(no_file, ew_order, EPIteNo, toc, ritz_ew(1), rsdl);
    
    % --- construct the initial search space
    %
    % Now, the first k columns of mtxV are the convergent eigenvectors.
    % The good initial vectors will be appended into mtxV from k+1-th column.
    %
    mtxV(:,1:k) = EV(:,1:k);
    j1 = 0;
    for i = 2:ewno 
          if ( flag_cn(i) == 0 )
             if ( j1+1 <= ew_number+RestartProjProbDim-1 ) 
                mtxV(:,k+j1+1) = init_Vec(:,i); 
             end
          end
          j1 = j1 + 1; 
    end

    fprintf('\n');
    
    %
    % Find the orthonormal basis of mtxV and reset it into init_vec
    %
    temp = mtxV;
    if  flag_no_ew <= Orth_Number
         [i, mtxV(:,1:i)]          = Gram_Schmidt(k+j1, mtxV(:, 1:k+j1));
         ProjProbDim               = min(i, ew_number+RestartProjProbDim-1) - k;
    elseif flag_no_ew > Orth_Number
    [i, mtxV(:,k-Orth_Number +1:k+j1)]          = Gram_Schmidt(Orth_Number+j1, mtxV(:, k-Orth_Number +1 : k+j1));
    ProjProbDim               = 3;
    end
%     %%modify by siu
%          [i, mtxV(:,1:i)]          = Gram_Schmidt(k+j1, mtxV(:, 1:k+j1));
%          ProjProbDim               = min(i, ew_number+RestartProjProbDim-1) - k;
%     [i, temp(:,1:i)]          = Gram_Schmidt1(k+j1, temp(:, 1:k+j1));
%     ProjProbDim               = min(i, ew_number+RestartProjProbDim-1) - k;
%     [k+j1 i ew_number+RestartProjProbDim-1 k ProjProbDim]
%     ProjProbDim               = 3;
    if ProjProbDim < 3 
        a = 0;
    end
    init_Vec(:,1:ProjProbDim) = mtxV(:,k+1:k+ProjProbDim);   %%% let init_Vec orthgal to EV
    ritz_ew                   = ritz_ew(2:ProjProbDim+1,1);

    ew_order                  = ew_order + j + 1; 
    ep_count = ew_order;
end

% =======================================
%
 function u = Amtrix(x, A)
   u = A * x;
 end
