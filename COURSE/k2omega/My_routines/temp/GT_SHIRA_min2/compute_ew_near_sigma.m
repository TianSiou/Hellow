function [ deflation_enforce, lambda_in, ew_new, ev_new, conv_radius, rsdl, no_deuplicate_ew, cong_flag, ew_unit_c, ev_on_c, ...
    ew_inside_c_new, ev_inside_c_new, no_ew_in_c, next_ref_ev_outside_c, no_ref_ev, idxx, no_unit_ew ] = ...
    compute_ew_near_sigma( dim_n, mtx_NME, sigma, deflation_enforce, rad, rad_old, ew_new, ...
    ev_new, eigenwanted, tolerance, ew_unit_c, no_ew_in_c, ew_inside_c_new, ev_on_c, ...
    ev_inside_c_new, next_ref_ev_outside_c, no_ref_ev, mtx_M, mtx_L, no_unit_ew )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

        deflate_tol        = 1.0e-9;
        maxit_GTSHIRA      = 30;
        
        W                  = mtx_NME.mtx_A - sigma*mtx_NME.mtx_Q + (sigma^2)*(mtx_NME.mtx_A.');
        LU_W.Perm_amd_vec  = amd(W);
        W_reorder          = W(LU_W.Perm_amd_vec, LU_W.Perm_amd_vec);
        LU_W.Perm_amd      = sparse(LU_W.Perm_amd_vec,1:dim_n,ones(dim_n,1));
        
        [LU_W.Low_L,LU_W.upper_U,LU_W.Perm_LU] = lu(W_reorder);
        
        LU_W.Perm_LU_Tran  = LU_W.Perm_LU.';
        LU_W.Low_L_Tran    = LU_W.Low_L.';
        LU_W.upper_U_Tran  = LU_W.upper_U.'; 
        LU_W.sigma         = sigma; 
        
        if ( abs(sigma) <= 3.0e-3 )
            LU_W.mtx_Qp  = W;
            LU_W.mtx_QpT = W.';
        end
            
        restart_info  = 'no';
        
            if (strcmp(deflation_enforce, 'yes') && rad > 6.0e-4 )
                
                % 
                % Find the eigenvalues to be in the convergent region
                %
%                 if ( ~isempty(target_ew_skew_Cayley) )
%                     max_rad           = max(abs(target_ew_skew_Cayley+1./target_ew_skew_Cayley-mu_0));
%                     tmp2              = ref_ew_inside_c + 1./ref_ew_inside_c; 
%                     idx_defl          = find(abs(tmp2 - mu_0) < 1.4 * max_rad);
%                     lambda_defl       = ref_ew_inside_c(idx_defl,1);
%                     ev_lambda_dfl     = ref_ev_inside_c(:,idx_defl);
%                     ev_lambda_inv_dfl = ref_ev_outside_c(:,idx_defl);
%                 
%                     lambda            = [ lambda_defl; ew_new(1:2:end,1) ];
%                     ev_lambda         = [ ev_lambda_dfl(1:dim_n,:)      ev_new(1:dim_n,1:2:end) ];
%                     ev_lambda_inv     = [ ev_lambda_inv_dfl(1:dim_n,:)  ev_new(1:dim_n,2:2:end) ];
%                 else
%                     lambda            = ew_new(1:2:end,1);
%                     ev_lambda         = ev_new(1:dim_n,1:2:end);
%                     ev_lambda_inv     = ev_new(1:dim_n,2:2:end);
%                 end

%                 idx_lambda    = find(abs(lambda) <= rad_old);

                if ( isempty(ew_new) )
                    idx_lambda = [];
                else
                    lambda        = ew_new(1:2:end,1);
                    idx_lambda    = find(abs(lambda) <= rad_old);
                end
                
                if ( isempty(idx_lambda) )
                    flag = 1;
                else
                    %lambda            = ew_new(1:2:end,1);
                    
                    no_defl_ew    = 20;
                    [~,idx_sort]  = sort(abs(lambda(idx_lambda) + 1./lambda(idx_lambda) - sigma - 1/sigma));
                    if ( length(idx_sort) >= no_defl_ew )
                        idx_lambda = idx_lambda(idx_sort(1:no_defl_ew));
                    end
                    
                    ev_lambda     = ev_new(1:dim_n,1:2:end);
                    ev_lambda_inv = ev_new(1:dim_n,2:2:end);
                    lambda        = lambda(idx_lambda,1);
                    ev_lambda     = ev_lambda(:,idx_lambda);
                    ev_lambda_inv = ev_lambda_inv(:,idx_lambda);
    
                    idx_real      = find(abs(imag(lambda)) <= 1.0e-10);
                    idx_cmp       = find(abs(imag(lambda)) > 1.0e-10 & imag(lambda) > 0);
                    lambda        = lambda([idx_real; idx_cmp]);
                    ev_lambda     = ev_lambda(:,[idx_real; idx_cmp]);
                    ev_lambda_inv = ev_lambda_inv(:,[idx_real; idx_cmp]);
                    fprintf('Deflation case with length = %2.0f \n', length(lambda));
                
                    if ( real(sigma) >= -0.01 )
                        [ mtx_deflate, flag ] = construct_deflated_mtx(lambda, ev_lambda, ev_lambda_inv, ...
                              mtx_NME, sigma, LU_W, deflate_tol); 
                    else
                        [ mtx_deflate, flag ] = construct_deflated_mtx_pone(lambda, ev_lambda, ev_lambda_inv, ...
                              mtx_NME, sigma, LU_W, deflate_tol);
                    end
                end
                
                if ( flag == 0 )
                    if ( real(sigma) >= -0.01 )
                        [ ew_new, ev_new, outer_it ] = Drive_GTSHIRA( dim_n, ...
                             @(x,transp_flag)Defl_mtx_vec_A( x, transp_flag, mtx_NME, mtx_deflate ), ...
                             @(x)Defl_mtx_vec_Q( x, mtx_NME, mtx_deflate ), mtx_NME, ...
                             @(x,transp_flag)solve_LS_PQEP(x,transp_flag,LU_W, mtx_deflate), ...
                             sigma, eigenwanted, tolerance, maxit_GTSHIRA, restart_info, [], [], [], [] );
                    else
                        
                        [ ew_new, ev_new, outer_it ] = Drive_GTSHIRA( dim_n, ...
                             @(x,transp_flag)Defl_mtx_vec_A( x, transp_flag, mtx_NME, mtx_deflate ), ...
                             @(x)Defl_mtx_vec_Q_pone( x, mtx_NME, mtx_deflate ), mtx_NME, ...
                             @(x,transp_flag)solve_LS_PQEP(x,transp_flag,LU_W, mtx_deflate), ...
                             sigma, eigenwanted, tolerance, maxit_GTSHIRA, restart_info, [], [], [], [] );
                    
                    end
                    
                    deflation_enforce = 'no';
                else
                    [ ew_new, ev_new, outer_it ] = Drive_GTSHIRA( dim_n, mtx_NME.mtx_A, mtx_NME.mtx_Q, mtx_NME, ...
                        @(x,transp_flag)solve_LS_PQEP(x,transp_flag,LU_W), sigma, eigenwanted, tolerance, maxit_GTSHIRA, ...
                        restart_info, [], [], [], [] );
                    
                    deflation_enforce = 'yes';
                end

            else
                [ ew_new, ev_new, outer_it ] = Drive_GTSHIRA( dim_n, mtx_NME.mtx_A, mtx_NME.mtx_Q, mtx_NME, ...
                    @(x,transp_flag)solve_LS_PQEP(x,transp_flag,LU_W), sigma, eigenwanted, tolerance, maxit_GTSHIRA, ...
                        restart_info, [], [], [], [] );
                
                if ( strcmp(deflation_enforce, 'no') )
                    deflation_enforce = 'yes';
                end

            end
            
            if ( ~isempty(ew_new) ) % GTSHIRA converges or not
                cong_flag     = 0;
                mu_0          = sigma + 1 / sigma;
                conv_radius   = max(abs(ew_new + 1./ew_new - mu_0));
                    
                %
                % Find the deuplicated eigenvalues
                %
                leng_ew       = length(ew_new);
                lambda_in     = ew_new(1:2:leng_ew);
                
                if ( ~isempty(find(abs(lambda_in) > 1+1.0e-7, 1)) )
                    fprintf('Error in choosing lambda_in \n'); 
                    return
                end
                
                ev_in_c       = ev_new(:,1:2:leng_ew);
                ev_out_c      = ev_new(:,2:2:leng_ew);
                idx_deupli    = zeros(leng_ew/2,1);
                idx_new_ui_ew = zeros(leng_ew/2,1);
             
                %
                % (i) Check any deuplicated or lost unimodular eigenvalues
                %
                tol        = 5.0e-6;
                idx_unit   = find(abs(lambda_in) <= 1+tol & abs(lambda_in) >= 1-tol);
                if ( ~isempty(idx_unit))
                    for i1 = 1:length(idx_unit)
                        if (min(abs(lambda_in(idx_unit(i1))-ew_unit_c(1:no_unit_ew,1))) < tol )
                            idx_deupli(idx_unit(i1),1) = 1; 
                        elseif (min(abs(lambda_in(idx_unit(i1))-ew_unit_c(1:no_unit_ew,2))) < tol )
                            idx_deupli(idx_unit(i1),1) = 1;
                        else
                            % Such unimodular eigenvalue is lost in the
                            % previous computations
                            idx_new_ui_ew(idx_unit(i1),1) = 1;
                        end
                    end
                end
            
                %
                % (ii) Check any deuplicated eigenvalues in the unit circle
                %
                idx_inside = setdiff(1:leng_ew/2, idx_unit); 
                if ( isempty(idx_inside) )
                    fprintf('all eigenvalues are unimodular eigenvalues\n');
                else
                    %
                    % Define which inside eigenvalues to be used checking the new computing eigenvalues
                    % being deuplicated or not
                    %
                    %interval = 1:no_ew_in_c;
%                     if ( jj < 3 )
%                         interval = 1:idx_ew_in_c(jj,1);
%                     else
%                         interval = idx_ew_in_c(jj-2,1):idx_ew_in_c(jj,1);
%                     end
                    interval = 1:no_ew_in_c;
                
                    for i1 = 1:length(idx_inside)
                        if (min(abs(lambda_in(idx_inside(i1))-ew_inside_c_new(interval,1)))/abs(lambda_in(idx_inside(i1))) < tol )
                            idx_deupli(idx_inside(i1),1) = 1;
                        end
                    end
                
                end
            
                idxx = find(idx_new_ui_ew == 1);
                if ( ~isempty(idxx) )
                    mm3                                      = length(idxx);  
                    
                    ew_unit_c(no_unit_ew+1:no_unit_ew+mm3,1) = lambda_in(idxx);
                    ew_unit_c(no_unit_ew+1:no_unit_ew+mm3,2) = ew_new(2*idxx);
                    ev_on_c(:,no_unit_ew+1:no_unit_ew+mm3,1) = ev_in_c(:,idxx);
                    ev_on_c(:,no_unit_ew+1:no_unit_ew+mm3,2) = ev_out_c(:,idxx);
                    no_unit_ew                               = no_unit_ew + mm3; 
                    
%                     fprintf('no_unit_ew = %3.0f in compute_ew_near_sigma \n', no_unit_ew);
%                     ew_unit_c(1:no_unit_ew,1)
                    
%                     ew_unit_c(mm+1:mm+length(idxx),1) = lambda_in(idxx); %[ew_unit_c(:,1); lambda_in(idxx)];
%                     ew_unit_c(mm+1:mm+length(idxx),2) = ew_new(idxx+1); %[ew_unit_c(:,2); ew_new(idxx+1)];
%                     ev_on_c(:,mm+1:mm+length(idxx),1) = ev_in_c(:,idxx);
%                     ev_on_c(:,mm+1:mm+length(idxx),2) = ev_out_c(:,idxx);
                end
            
                idxx                   = find(idx_deupli == 0); 
                mm2                    = length(idxx);
                percentage             = (leng_ew/2-mm2)/(leng_ew/2);
                no_deuplicate_ew       = [ leng_ew/2-mm2  (leng_ew/2-mm2)/(leng_ew/2)];
                fprintf('Re_GTSHIRA, percentage of deuplicated eigenvalues = %11.4e \n', percentage)
                %if ( mm2 ~= leng_ew/2 )
                    %
                    idx1              = find(imag(lambda_in(idxx)) > -1.0e-10);
                    mm2               = length(idx1);
                    idxx              = idxx(idx1);
                    %
                    idx2              = zeros(mm2*2,1); 
                    idx2(1:2:2*mm2,1) = 2*idxx-1;
                    idx2(2:2:2*mm2,1) = 2*idxx;
                    ew_new            = ew_new(idx2);
                    ev_new            = ev_new(:,idx2);
                %end
            
%                 ew_inside_c_new(no_ew_in_c+1:no_ew_in_c+mm2,1)     = lambda_in(idxx);
%                 ev_inside_c_new(:,no_ew_in_c+1:no_ew_in_c+mm2)     = ev_in_c(:,idxx);
%                 no_ew_in_c                                         = no_ew_in_c + mm2;  
%                 next_ref_ev_outside_c(:,no_ref_ev+1:no_ref_ev+mm2) = ev_out_c(:,idxx);
%                 no_ref_ev                                          = no_ref_ev + mm2;
                
                rsdl = zeros(length(ew_new)/2,1);
                for ki = 1:2:length(ew_new)
                    if abs(imag(ew_new(ki)))<1e-10
                        ew_new(ki, 1) = real(ew_new(ki, 1));
                        ev_new(:, ki) = real(ev_new(:, ki));
                    end
                    ev_new(:,ki)     = ev_new(:,ki) / norm(ev_new(:,ki));
                    rsdl_vec         = mtx_M * ev_new(:,ki) - ew_new(ki) * (mtx_L * ev_new(:,ki)); 
                    rsdl((ki+1)/2,1) = norm(rsdl_vec);
                    fprintf('rsdl(%2.0f,%24.16e+(%24.16e)1i) = %10.4e \n', (ki+1)/2, real(ew_new(ki)), ...
                    imag(ew_new(ki)), rsdl((ki+1)/2,1)); 
                end
                
                idx_rsdl          = find( rsdl <= 1.0e-9 );
                leng_idx_rsdl     = length(idx_rsdl);
                if ( leng_idx_rsdl == length(ew_new)/2 )
                    ew_inside_c_new(no_ew_in_c+1:no_ew_in_c+mm2,1)     = lambda_in(idxx);
                    ev_inside_c_new(:,no_ew_in_c+1:no_ew_in_c+mm2)     = ev_in_c(:,idxx);
                    no_ew_in_c                                         = no_ew_in_c + mm2;  
                    next_ref_ev_outside_c(:,no_ref_ev+1:no_ref_ev+mm2) = ev_out_c(:,idxx);
                    no_ref_ev                                          = no_ref_ev + mm2;
                elseif ( leng_idx_rsdl > 0 )
                    ew_inside_c_new(no_ew_in_c+1:no_ew_in_c+leng_idx_rsdl,1)     = lambda_in(idxx(idx_rsdl));
                    ev_inside_c_new(:,no_ew_in_c+1:no_ew_in_c+leng_idx_rsdl)     = ev_in_c(:,idxx(idx_rsdl));
                    no_ew_in_c                                                   = no_ew_in_c + leng_idx_rsdl;  
                    next_ref_ev_outside_c(:,no_ref_ev+1:no_ref_ev+leng_idx_rsdl) = ev_out_c(:,idxx(idx_rsdl));
                    no_ref_ev                                                    = no_ref_ev + leng_idx_rsdl;
                    
                    % == add 2016/12/26
                    lambda_in         = lambda_in(idxx(idx_rsdl)); 
                    
                    mm2               = length(idx_rsdl);
                    idx2              = zeros(mm2*2,1); 
                    idx2(1:2:2*mm2,1) = 2*idx_rsdl-1;
                    idx2(2:2:2*mm2,1) = 2*idx_rsdl;
                    ew_new            = ew_new(idx2);
                    ev_new            = ev_new(:,idx2);
                    idxx              = 1:mm2; 
                else
                    cong_flag    = 1;
                end
                
            else
                cong_flag        = 1;
                conv_radius      = 0; 
                rsdl             = []; 
                no_deuplicate_ew = [];
                lambda_in        = [];
                idxx             = [];
            end
                
end

