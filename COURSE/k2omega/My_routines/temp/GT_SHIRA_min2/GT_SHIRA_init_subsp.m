function [ev, ew, Hess_H, UpTr_R, inv_Z, outer_it] = GT_SHIRA_init_subsp(n, mtx_A, mtx_Q, solve_LS_PQEP, ... 
        xini, eigswanted, nstp, itmax, tolerance, shift, dynamical_change_eigswanted, init_R, init_H, init_Y, init_Z) 
%
%    This program is a solver for finding the eigenvalues with the largest
%  modular for the eigenproblem K x = mu N x where K, N are
%  T-skew-Hamiltonian. The program is applied 
%  by (T-skew-Hamiltonian) structure-preserving Arnoldi algorithm.
%
%  Input:
%        Matrices K, N with dim = 2n.
%          xini = the initial vector.
%          nini = the choice parameter for determining the initial vector.
%        maxite = upper bound of Lanczos iteration.
%            no = total number of convergent eigenvales you want.
%
%  Output:
%         ev  = convergence right eigenvectors
%         ew  = convergence eigenvalues
%         err = error message
%
%  Notice:
%         'err' is an important parameter for this program because
%         the algorithm needs restart when an invariant subspace ocurrs.
%

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% call skew-Hamiltonian Arnoldi iteration
% 

nstp_org = nstp;
n2       = 2 * n; 
%nmax     = 5*nstp;  %5*nstp;  %4*nstp + 5;
if eigswanted < 20
    nmax     = 5 * nstp;
elseif eigswanted < 100
    nmax     = 3 * nstp;
else
    nmax     = 2 * nstp;
end
fprintf('number of wanted eigenvalues = %3.0f, restart number = %2.0f \n',eigswanted, nmax);
R        = zeros(nmax);
H        = zeros(nmax+1,nmax);
Y        = zeros(n2,nmax+1);
Z        = zeros(n2,nmax); 
% ew       = 0;
% ev       = 0;

%if eigswanted > nstp - 1
%  eigswanted = nstp - 1;
%end

%
% Generate starting vector for Arnoldi process.
%  
Y(:,1) = xini; 
Y(:,1) = Y(:,1) / norm(Y(:,1));

% Do nstp Arnoldi steps.

loopstart    = 1; loopend = nstp; 
Y(:,1:loopend+1)         = init_Y;
Z(:,1:loopend)           = init_Z;
H(1:loopend+1,1:loopend) = init_H;
R(1:loopend,1:loopend)   = init_R;

% [Y(:,1:loopend+1), Z(:,1:loopend), H(1:loopend+1,1:loopend), R(1:loopend,1:loopend)] =  ... 
%      GTArnoldi_step(n, mtx_A, mtx_Q, solve_LS_PQEP, Y(:,1:loopstart), Z(:,1:loopstart-1), H(1:loopstart,1:loopstart-1), ...
%      R(1:loopstart-1,1:loopstart-1), shift, loopstart, loopend);

% KZ    = mtx_K_prod_vec( mtx_G, mtx_M2, mtx_F, mtx_M1, Z(:,1:loopend), shift );
% Err_K = KZ - Y(:,1:loopend) * H(1:loopend,1:loopend) - H(1+loopend,loopend) * Y(:,loopend+1) * [ zeros(1,loopend-1) 1];
% fprintf('Err_invK = %14.4e \n', norm(Err_K));
% NZ    = mtx_N_prod_vec( mtx_G, mtx_M2, mtx_F, mtx_M1, Z(:,1:loopend), shift );
% Err_N = NZ - Y(:,1:loopend) * R(1:loopend,1:loopend);
% fprintf('Err_invN = %14.4e \n', norm(Err_N));
% if ( flag ~= 0 ) 
%     return
% end

flag_change_no_ew = 'yes';
stop_change_no_ew = 'no';
converged         = 0;
no_change_no_ew   = 0;

conv_flag         = 0;
j                 = 1;
Error_check       = 0;
no_conv_ew        = 0;
count_no_conv_ew  = 0;

while ( j <= itmax && conv_flag == 0 )
%for j = 1:itmax
   fprintf('%3.0f-th iteration of GTSHIRA with %2.0f convergent eigenvalues \n',j,converged);
   % Check for convergence.
   for k = converged+1:nstp-1
      if abs(H(k+1,k)) < tolerance*(abs(H(k,k) + H(k+1,k+1)))
         converged           = k; 
         H(k+1,k)            = 0;
      end
   end
   if converged < eigswanted
       [REV,REW]       = eig(H(1:nstp,1:nstp), R(1:nstp,1:nstp));
       residual        = abs(H(nstp+1,nstp)) * abs(REV(nstp,:));
       idx_conv        = find(residual <= tolerance);
       prev_no_conv_ew = no_conv_ew;
       no_conv_ew      = length(idx_conv);
       if ( prev_no_conv_ew == no_conv_ew && no_conv_ew > 0 )
           count_no_conv_ew = count_no_conv_ew + 1;
       else
           count_no_conv_ew = 0;
       end
       fprintf('nstp = %2.0f, number of convergent eigenvalues = %2.0f \n',nstp,no_conv_ew);
       
       if (no_conv_ew == nstp)
           conv_flag = 1;
           converged = nstp;
       else
           if ( strcmp(dynamical_change_eigswanted,'yes') || strcmp(dynamical_change_eigswanted,'YES') )
               
               if ( j >= min(80, itmax-5) )
%                    nstp                        = max(converged,no_conv_ew);
%                    eigswanted                  = min(eigswanted, nstp);
                   dynamical_change_eigswanted = 'no';
               elseif ( strcmp(stop_change_no_ew,'no') && j > 15 && strcmp(flag_change_no_ew,'yes') ) 
                   if ( nstp > eigswanted )
                       ew_sort           = diag(REW);
                       [~,idx_sort]      = sort(-abs(ew_sort));
                       ew_sort           = ew_sort(idx_sort);
                       diff_ew           = abs(ew_sort(1:end-1) - ew_sort(2:end));
                       [~,idx_max_gap]   = max(diff_ew(eigswanted:end));
                       no_dyn_ew         = eigswanted - 1 + idx_max_gap;
                       nstp              = max(no_dyn_ew, no_conv_ew);
                       no_change_no_ew   = no_change_no_ew + 1;
                       fprintf('new nstp = %2.0f \n',nstp);
                   end
                   flag_change_no_ew = 'no';
               
               elseif ( strcmp(stop_change_no_ew,'no') && j >= 20*no_change_no_ew && strcmp(flag_change_no_ew,'no') )
                   flag_change_no_ew = 'yes';
               
               elseif ( count_no_conv_ew >= 10 )
                   nstp       = no_conv_ew;
                   eigswanted = min(eigswanted, nstp);
               end
               
           end

           % Do another nstp+1 Arnoldi steps.

           loopstart  = nstp_org+1; loopend = nmax;
           [Y(:,1:loopend+1), Z(:,1:loopend), H(1:loopend+1,1:loopend), R(1:loopend,1:loopend)] =  ... 
               GTArnoldi_step(n, mtx_A, mtx_Q, solve_LS_PQEP, Y(:,1:loopstart), Z(:,1:loopstart-1), H(1:loopstart,1:loopstart-1), ...
               R(1:loopstart-1,1:loopstart-1), shift, loopstart, loopend);
   
           nstp_org = nstp;
   
           if ( Error_check == 1 )
               fprintf('Error check for Krylov subspace before restarting \n')
               KZ    = mtx_K_prod_vec( mtx_A, Z(:,1:loopend), shift );
               Err_K = KZ - Y(:,1:loopend) * H(1:loopend,1:loopend) - ...
                     H(1+loopend,loopend) * Y(:,loopend+1) * [ zeros(1,loopend-1) 1];
               fprintf('Error of Arnoldi decomp with K = %11.4e\n', norm(Err_K));
               NZ    = mtx_N_prod_vec( mtx_A, mtx_Q, Z(:,1:loopend), shift );
               Err_N = NZ - Y(:,1:loopend) * R(1:loopend,1:loopend);
               fprintf('Error of Arnoldi decomp with N = %11.4e\n',norm(Err_N,inf))
           end

           % Retain the first nstp columns.

           Restart_method = 'Schur'; %'Implicit'; %'Schur';
           switch Restart_method
               case 'Implicit'
                   tmp_H                = H(loopend+1,loopend);
                   test_conv            = 0; %converged
%           [H, R, Y, Z, unit_e] = Refine_shift_Implicit_restart(test_conv, nmax, nstp, H, R, Y, Z); 
                   [H, R, Y, Z, unit_e] = Implicit_restart(test_conv, nmax, nstp, H, R, Y, Z);
                   H(nstp+1,nstp)       = tmp_H * unit_e(nstp);
               case 'Schur'
                   tmp_H                   = H(loopend+1,loopend);
                   vec_y                   = Y(:,loopend+1);
                      
                   [new_H, new_R, Y_new, Z_new, alpha, new_nstp] = Schur_restart(nmax, nstp_org, H, R, Y, Z);
           
                   R(1:new_nstp,1:new_nstp) = new_R;
                   H(1:new_nstp,1:new_nstp) = new_H;
                   H(1+new_nstp,1:new_nstp) = [ zeros(1,new_nstp-1) alpha*tmp_H ];
                   Y(:,1:new_nstp)          = Y_new;
                   Y(:,1+new_nstp)          = vec_y;
                   Z(:,1:new_nstp)          = Z_new;
                   nstp                     = new_nstp;
           end

           if ( Error_check == 1 )
               fprintf('Error check after restarting \n')
               KZ    = mtx_K_prod_vec( mtx_A, Z(:,1:nstp), shift );
               Err_K = KZ - Y(:,1:nstp) * H(1:nstp,1:nstp) - H(1+nstp,nstp) * Y(:,nstp+1) * [ zeros(1,nstp-1) 1];
               fprintf('Error of Arnoldi decomp with K = %11.4e\n', norm(Err_K));
               NZ    = mtx_N_prod_vec( mtx_A, mtx_Q, Z(:,1:nstp), shift );
               Err_N = NZ - Y(:,1:nstp) * R(1:nstp,1:nstp);
               fprintf('Error of Arnoldi decomp with N = %11.4e\n',norm(Err_N,inf))
           end

           j = j + 1;
       end
   else
       conv_flag = 1;
       
   end
   
end

   if ( converged >= eigswanted || conv_flag == 1 ) % Compute eigenvalues and stop.
%       [ew, ev]     = irafinish(converged, H, R, Z);
      [REV,REW]    = eig(H(1:converged,1:converged), R(1:converged,1:converged));
      ew           = diag(REW);
      refine_flag  = 0;
      if ( refine_flag == 0 )
          ev = Z(:,1:converged) * REV;
      else
%           ev = refinement(H(1:converged,1:converged), R(1:converged,1:converged), ew, Z(:,1:converged));
          [ew, ev] = refinement(H(1:converged+1,1:converged), R(1:converged+1,1:converged), ew, Z(:,1:converged));
      end
      inv_Z        = Z(:,1:converged);
      Hess_H       = H(1:converged,1:converged);
      UpTr_R       = R(1:converged,1:converged);
      outer_it     = j; 
   else
       [REV,REW]   = eig(H(1:nstp,1:nstp), R(1:nstp,1:nstp));
       residual    = abs(H(nstp+1,nstp)) * abs(REV(nstp,:));
       idx_conv    = find(residual <= tolerance);
       if ( ~isempty(idx_conv) )
           ew          = diag(REW);
           ew          = ew(idx_conv);
           ev          = Z(:,1:nstp) * REV(:,idx_conv);
       else
           ew = [];
           ev = [];
       end
       inv_Z       = [];
       Hess_H      = [];
       UpTr_R      = [];
       outer_it    = j;
   end
   
%
% if ( flag == 0 )
%     flag = -1;
% end
% disp('Not Done!')


