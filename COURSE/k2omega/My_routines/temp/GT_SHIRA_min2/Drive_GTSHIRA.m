function [ ew_new, ev_new, outer_it, convg_lambda, counting_ew, no_deuplicate_ew ] = Drive_GTSHIRA( n, mtx_A, ...
    mtx_Q, mtx_NME, solve_LS_PQEP, shift, eigenwanted, tolerance, maxit, ...
    restart_info, R, H, Y, Z, convg_lambda, counting_ew)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
% OUTPUT
%
%   ew_new      : one dimensional array
%                 ew_new(1:2:end) save the eigenvalues lambda with abs(lambda) <= 1
%                 ew_new(2:2:end) save the eigenvalues lambda with abs(lambda) >= 1
%                 ew_new(2*i-1) * ew_new(2*i) = 1
%   ev_new      : two-dimensional array
%                 ev_new(:,1:2:end) save the eigenvectors corresponding to ew_new(1:2:end)
%                 ev_new(:,2:2:end) save the eigenvectors corresponding to ew_new(2:2:end)
%   ew_unit_c   : two-dimensional array
%                 ew_unit_c(:,1) and ew_unit_c(:,2) save the unimodular eigenvalues
%                 ew_unit_c(i,1) * ew_unit_c(i,2) = 1
%   ew_inside_c : one-dimensional array
%                 save the eigenvalues which are inside the circle
%   convg_lambda: two-dimensional array
%                 the mm-th column is used to save the convergent eigenvalues in the current step
%                 the (mm-1)-th column is used to save the convergent eigenvalues in the previous step
%                 and so on
%   counting_ew : one dimensional array (corresponding to convg_lambda)
%                 the mm-th element is used to save the number of the convergent eigenvalues in the current step
%                 the (mm-1)-th element is used to save the number of the convergent eigenvalues in the previous step
%                 and so on

    nstp  =   eigenwanted; %+ 10;
%     maxit = 100;
    
    s = RandStream('mt19937ar','Seed',0);
    RandStream.setGlobalStream(s);
    xini  = randn(2*n,1) + 1i * randn(2*n,1);
    xini  = xini / norm(xini);
    
    dynamical_change_eigswanted = 'yes';

%     W                  = mtx_A - shift*mtx_Q + (shift^2)*(mtx_A.');
%     LU_W.Perm_amd_vec  = amd(W);
%     W_reorder          = W(LU_W.Perm_amd_vec, LU_W.Perm_amd_vec);
%     LU_W.Perm_amd      = sparse(LU_W.Perm_amd_vec,1:n,ones(n,1));
% 
%     [LU_W.L,LU_W.U,LU_W.Perm_LU] = lu(W_reorder);
    
%     [REV, REW, H, R, Z, converged] = GT_SHIRA(n, mtx_A, mtx_Q, LU_W, xini, eigenwanted, ...

    if ( strcmp(restart_info,'no') )
        [REV, REW, ~, ~, ~, outer_it] = GT_SHIRA(n, mtx_A, mtx_Q, solve_LS_PQEP, xini, eigenwanted, ...
             nstp, maxit, tolerance, shift, dynamical_change_eigswanted); 
    else
%         [REV, REW, ~, ~, ~, outer_it] = GT_SHIRA(n, mtx_A, mtx_Q, solve_LS_PQEP, xini, eigenwanted, ...
%              nstp, maxit, tolerance, shift, dynamical_change_eigswanted); 
        nstp = size(R,1);
        [REV, REW, ~, ~, ~, outer_it] = GT_SHIRA_init_subsp(n, mtx_A, mtx_Q, solve_LS_PQEP, xini, eigenwanted, ...
             nstp, maxit, tolerance, shift, dynamical_change_eigswanted, R, H, Y, Z);
    end
    
    if ( isempty(REW) )
        ew_new           = []; 
        ev_new           = []; 
        outer_it         = 0; 
        convg_lambda     = []; 
        counting_ew      = 0; 
        no_deuplicate_ew = 0;
    else
        ll         = size(REW  ,    1);  

        %
        % Compute the eigenvalues lambda_vec and lambda_vec.^(-1) of the symplectic pencil (M, L)
        %
        lambda_vec = shift + (1.0/shift) + 1./REW;
        tmp        = lambda_vec.^2 - 4;
        lambda_vec = (lambda_vec + tmp.^(0.5))/2;
    
        lambda_tmp = lambda_vec;
        tol        = 1.0e-7;
        for jj = 1:length(lambda_vec) 
            lambda = lambda_vec(jj,1);
            if ( abs(lambda) > 1+tol )
                lambda_tmp(jj, 1) = 1/lambda; 
            elseif ( abs(lambda) < 1-tol )
                lambda_tmp(jj, 1) = lambda  ; 
            else
                if ( imag(lambda) > 0 )
                    lambda_tmp(jj, 1) = lambda  ; 
                else
                    lambda_tmp(jj, 1) = 1/lambda; 
                end
            end 
        end
    
        [~,idx_sort]          = sort(real(lambda_tmp));
        lambda_tmp            = lambda_tmp(idx_sort);
    
        lambda_vec            = lambda_vec(idx_sort);
        REV                   = REV(:,idx_sort);
    
        if ( nargin > 14 && nargout > 3 )
            
%             figure(10)
%             theta = 0:0.001:2*pi;
%             x     = cos(theta)+1i*sin(theta);
% 
%             plot(x,'r-')
%             axis([-1.1, 1.1, -1.1, 1.1]);
%             hold on
%             load ../v1_0919/EW_data_n2p0
%             plot(ew,'bo','LineWidth',1);
%             plot(lambda_tmp,'g+','LineWidth',2);
%             hold off 
            
            deuplicate_lambda = zeros(ll,1);
            mm                = size(convg_lambda,2);
            for ii = mm:-1:1
                if (counting_ew(ii) > 0)
                    idx = find(real(lambda_tmp) > real(convg_lambda(1,ii))-1.0e-8); 
                    for i2 = 1:length(idx)
                        if (min(abs(lambda_tmp(idx(i2))-convg_lambda(:,ii))) < 1.0e-7 )
                            deuplicate_lambda(idx(i2),1) = 1;
                            fprintf('deuplication ew(%2.0f) = %11.4e+(%11.4e)i \n', ii, ...
                                real(lambda_tmp(idx(i2))), imag(lambda_tmp(idx(i2))))
                        end
                    end
                end
            end
    
            idx              = find(deuplicate_lambda == 0);
            no_deuplicate_ew = [ll - length(idx) ll];
    
            lambda_vec       = lambda_vec(idx);
            REV              = REV(:,idx);
            lambda_tmp       = lambda_tmp(idx);
    
            convg_lambda(:,1:mm-1) = convg_lambda(:,2:mm);
            i2                     = length(lambda_tmp);
            ij                     = size(convg_lambda,1) - i2;
            if ( ij > 0 )
                convg_lambda(:,mm)    = [ lambda_tmp; zeros(ij,1) ];
            else
                convg_lambda(1:i2,mm) = lambda_tmp;
            end
    
            counting_ew(1:mm-1,1) = counting_ew(2:mm,1);
            counting_ew(mm,1)     = i2;
        
        end
    
        ll          = length(lambda_vec);
        ew_new      = zeros(2*ll,    1);
        ev_new      = zeros(2*n , 2*ll);
%     ew_unit_c   = zeros(2*ll,    2);
%     ew_inside_c = zeros(2*ll,    1);
%     no_uc       = 0;
%     no_ic       = 0;
    
        for jj = 1:length(lambda_vec) 
            lambda = lambda_vec(jj,1);
            if ( abs(lambda) > 1+tol )
                ew_new(2*jj-1, 1) = 1/lambda;
                ew_new(2*jj  , 1) = lambda  ;
            elseif ( abs(lambda) < 1-tol )
                ew_new(2*jj-1, 1) = lambda  ;
                ew_new(2*jj  , 1) = 1/lambda; 
%         elseif ( abs(lambda) <= 1 )
%             if ( imag(lambda) > 0 )
%                 ew_new(2*jj-1, 1) = lambda  ;
%                 ew_new(2*jj  , 1) = 1/lambda;
%             else
%                 ew_new(2*jj-1, 1) = 1/lambda;
%                 ew_new(2*jj  , 1) = lambda  ;
%             end
            else
                if ( imag(lambda) > 0 )
                    ew_new(2*jj-1, 1) = lambda  ;
                    ew_new(2*jj  , 1) = 1/lambda;
                else
                    ew_new(2*jj-1, 1) = 1/lambda;
                    ew_new(2*jj  , 1) = lambda  ;
                end
            end 
            ev_new(1:n    , 2*jj-1) = ew_new(2*jj, 1)*REV(1:n, jj) - REV(n+1:2*n, jj);
        
            ttQ = mtx_NME.mtx_Q * REV(1:n, jj);
        
            ttA  = mtx_NME.mtx_A * REV(1:n, jj);
            ttAT = mtx_NME.mtx_A.' * REV(n+1:2*n, jj);
 
            ev_new(n+1:2*n, 2*jj-1) = -ttA + ew_new(2*jj, 1) * ttQ - ew_new(2*jj, 1)*ttAT;
% %         ev_new(n+1:2*n, 2*jj-1) = -mtx_A * REV(1:n, jj) + ew_new(2*jj, 1) * (mtx_Q * REV(1:n, jj)) ...
% %             - ew_new(2*jj, 1)*(mtx_A' * REV(n+1:2*n, jj));

            ev_new(1:n    , 2*jj  ) = ew_new(2*jj-1, 1)*REV(1:n, jj) - REV(n+1:2*n, jj);
            ev_new(n+1:2*n, 2*jj  ) = -ttA + ew_new(2*jj-1, 1) * ttQ - ew_new(2*jj-1, 1)*ttAT;
% 
% %         ev_new(n+1:2*n, 2*jj  ) = -mtx_A * REV(1:n, jj) + ew_new(2*jj-1, 1) * (mtx_Q * REV(1:n, jj)) ...
% %             - ew_new(2*jj-1, 1)*(mtx_A' * REV(n+1:2*n, jj));
        
%         if ( abs(lambda) > 1+tol )
%             no_ic                = no_ic + 1;
%             ew_inside_c(no_ic,1) = 1 / lambda;
%         elseif ( abs(lambda) < 1-tol )
%             no_ic                = no_ic + 1;
%             ew_inside_c(no_ic,1) = lambda;
%         elseif ( abs(lambda) <= 1 )
%             no_uc              = no_uc + 1;
%             ew_unit_c(no_uc,1) = ew_new(2*jj-1, 1); %lambda;
%             ew_unit_c(no_uc,2) = ew_new(2*jj  , 1); %1 / lambda;
%             
%             if ( imag(ew_unit_c(no_uc,1)) >= 0 )
%                 if isa(mtx_A,'function_handle')
%                     ttAT = mtx_A(ev_new(1:n,2*jj-1), 'transp');
%                 else
%                     ttAT = mtx_A.' * ev_new(1:n,2*jj-1);
%                 end
%                 tt = ev_new(1:n,2*jj-1)' * ttAT;
%                 tt = ew_unit_c(no_uc,1) * tt;
%             else
%                 if isa(mtx_A,'function_handle')
%                     ttAT = mtx_A(ev_new(1:n,2*jj), 'transp');
%                 else
%                     ttAT = mtx_A.' * ev_new(1:n,2*jj);
%                 end
%                 tt = ev_new(1:n,2*jj)' * ttAT;
%                 tt = ew_unit_c(no_uc,2) * tt;
%             end
%             fprintf('lambda_xHATx = %10.4e+(%10.4e)*1i\n',real(tt),imag(tt));
%         else
%             no_uc              = no_uc + 1;
%             ew_unit_c(no_uc,1) = ew_new(2*jj-1, 1); %1 / lambda;
%             ew_unit_c(no_uc,2) = ew_new(2*jj  , 1); %lambda;
%             
%             if ( imag(ew_unit_c(no_uc,1)) >= 0 )
%                 if isa(mtx_A,'function_handle')
%                     ttAT = mtx_A(ev_new(1:n,2*jj-1), 'transp');
%                 else
%                     ttAT = mtx_A.' * ev_new(1:n,2*jj-1);
%                 end
%                 tt = ev_new(1:n,2*jj-1)' * ttAT;
%                 tt = ew_unit_c(no_uc,1) * tt;
%             else
%                 if isa(mtx_A,'function_handle')
%                     ttAT = mtx_A(ev_new(1:n,2*jj), 'transp');
%                 else
%                     ttAT = mtx_A.' * ev_new(1:n,2*jj);
%                 end
%                 tt = ev_new(1:n,2*jj)' * ttAT;
%                 tt = ew_unit_c(no_uc,2) * tt;
%             end
%             fprintf('lambda_xHATx = %10.4e+(%10.4e)*1i\n',real(tt),imag(tt));
%         end
        
        end
    end

%     if ( no_ic == 0 )
%         ew_inside_c = [];
%     else
%         ew_inside_c = ew_inside_c(1:no_ic,1);
%     end
%     
%     if ( no_uc == 0 )
%         ew_unit_c = [];
%     else
%         ew_unit_c = ew_unit_c(1:no_uc,:);
%     end
    
%     [~, ndx] = sort(abs(ew_new), 'descend');
%     ew_new = ew_new(ndx, 1);
%     ev_new = ev_new(:, ndx);

end

