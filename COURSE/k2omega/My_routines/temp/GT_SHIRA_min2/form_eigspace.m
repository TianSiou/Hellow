function [ X_sol, RelErr ] = form_eigspace( mtx_NME, mtx_M, mtx_L, ew_inside_c, ev_inside_c, ev_all_outside_c, ew_unit_c, ev_on_c )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% %n           = size(ev_inside_c,1);
% no_ew       = length(ew_inside_c);
% idx_unit_ew = find(abs(abs(ew_inside_c)-1) < 1.0e-10);
% idx_cmp_ew  = find(abs(imag(ew_inside_c)) > 1.0e-9 & abs(abs(ew_inside_c)-1) >= 1.0e-10);
% idx_real_ew = setdiff(1:no_ew, [idx_cmp_ew; idx_unit_ew]);
% eig_space   = [ real(ev_inside_c(:,idx_real_ew)), real(ev_inside_c(:,idx_cmp_ew)), imag(ev_inside_c(:,idx_cmp_ew)), ...
%     real(ev_inside_c(:,idx_unit_ew)), imag(ev_inside_c(:,idx_unit_ew)), real(ev_all_outside_c(:,idx_real_ew)), ...
%     real(ev_all_outside_c(:,idx_cmp_ew)), imag(ev_all_outside_c(:,idx_cmp_ew)) ];
% 
% % Err = eig_space' * [ eig_space(n/2+1:n,:); - eig_space(1:n/2,:) ];
% % 
% % nrm_1 = zeros(no_ew,1);
% % for ii = 1:no_ew
% %     nrm_1(ii,1) = norm(Err(:,ii),inf);
% % end

        s = RandStream('mt19937ar','Seed',0);
        RandStream.setGlobalStream(s);

        RelErr       = 1;
        tolerance    = 1.0e-12;
        dim_large    = size(ev_inside_c,1);
        dim_n        = dim_large / 2;
        no_ew        = length(ew_inside_c);
        tt           = [ew_unit_c(:,1); ew_unit_c(:,2)];
        %idx_unit_ew  = find( imag(tt) >= -1.0e-7 );
        idx_unit_ew2 = find(abs(abs(ew_inside_c)-1) < 1.0e-10);
        idx_cmp_ew   = find(abs(imag(ew_inside_c)) > 1.0e-9 & abs(abs(ew_inside_c)-1) >= 1.0e-10);
        idx_real_ew  = setdiff(1:no_ew, [idx_cmp_ew; idx_unit_ew2]);
        mm_real      = length(idx_real_ew);
        mm_cmp       = length(idx_cmp_ew);
        mm_unit_ew   = size(ew_unit_c,1); %length(idx_unit_ew);
        mm           = mm_real + 2 * mm_cmp;
        U_in         = zeros(dim_large,mm);
        U_out        = zeros(dim_large,mm);
        no_cong_ew   = mm_real + 2 * mm_cmp + mm_unit_ew;
        no_ew_in_c   = mm_unit_ew;
                
        U_in (:,        1:mm_real) = real(ev_inside_c(:,idx_real_ew));
        U_out(:,        1:mm_real) = real(ev_all_outside_c(:,idx_real_ew));
        for ii = 1:mm_real
            U_in(:,ii)  = U_in(:,ii) / norm(U_in(:,ii));
            U_out(:,ii) = U_out(:,ii) / (U_in(:,ii)' * [U_out(dim_n+1:dim_large,ii); -U_out(1:dim_n,ii)]);
%             rsdl1(ii) = norm(mtx_M * U_in(:,ii) - ew_inside_c(idx_real_ew(ii)) * (mtx_L * U_in(:,ii)))/abs(ew_inside_c(idx_real_ew(ii)));
%             rsdl2(ii) = norm(mtx_M * U_out(:,ii) - 1/ew_inside_c(idx_real_ew(ii)) * (mtx_L * U_out(:,ii)))/norm(U_out(:,ii))*abs(ew_inside_c(idx_real_ew(ii)));
        end
        kk                         = mm_real;
        
        for ii = 1:mm_cmp
            [ U_in(:,kk+1:kk+2), U_out(:,kk+1:kk+2) ] = J_orthogonal( [real(ev_inside_c(:,idx_cmp_ew(ii))) ...
                imag(ev_inside_c(:,idx_cmp_ew(ii)))], [real(ev_all_outside_c(:,idx_cmp_ew(ii))) ...
                imag(ev_all_outside_c(:,idx_cmp_ew(ii)))], dim_n );
            
%             rsdl1(ii) = norm(mtx_M * ev_inside_c(:,idx_cmp_ew(ii)) - ew_inside_c(idx_cmp_ew(ii)) * (mtx_L * ev_inside_c(:,idx_cmp_ew(ii))))/norm(ev_inside_c(:,idx_cmp_ew(ii)))/abs(ew_inside_c(idx_cmp_ew(ii)));
%             rsdl2(ii) = norm(mtx_M * ev_all_outside_c(:,idx_cmp_ew(ii)) - 1/ew_inside_c(idx_cmp_ew(ii)) * (mtx_L * ev_all_outside_c(:,idx_cmp_ew(ii))))/norm(ev_all_outside_c(:,idx_cmp_ew(ii)))*abs(ew_inside_c(idx_cmp_ew(ii)));
            
            kk                                        = kk + 2;
        end

    find_all_eigenvalue = 'yes';
    
    if ( no_cong_ew < dim_n )   
        %
        % Find lossing eigenpairs
        %
        find_all_eigenvalue = 'no';
        
        U_on_u = [real(ev_on_c(:,:,1)), imag(ev_on_c(:,:,1))];
        %U_on_u = [real(ev_inside_c(:,idx_unit_ew)), imag(ev_inside_c(:,idx_unit_ew))];
        for ii = 1:mm_unit_ew
            U_on_u(:,ii)            = U_on_u(:,ii) / norm(U_on_u(:,ii));
            U_on_u(:,ii+mm_unit_ew) = U_on_u(:,ii+mm_unit_ew) / (U_on_u(:,ii)'*[U_on_u(dim_n+1:dim_large,ii+mm_unit_ew); -U_on_u(1:dim_n,ii+mm_unit_ew)]);
            
%             rsdl1(ii) = norm(mtx_M * ev_on_c(:,ii,1) - ew_unit_c(ii,1) * (mtx_L * ev_on_c(:,ii,1)))/abs(ew_unit_c(ii,1))/norm(ev_on_c(:,ii,1));
%             rsdl2(ii) = norm(mtx_M * ev_on_c(:,ii,2) - ew_unit_c(ii,2) * (mtx_L * ev_on_c(:,ii,2)))/norm(ev_on_c(:,ii,2))/abs(ew_unit_c(ii,2));            
        end
        
        U_eig_space = [ U_in U_on_u(:,1:mm_unit_ew) U_out  U_on_u(:,1+mm_unit_ew:2*mm_unit_ew) ];
        
        checking_error = 'no';
        if ( strcmp(checking_error,'yes') )
            Err = U_in' * [ U_in(dim_n+1:dim_large,:); - U_in(1:dim_n,:) ];

            nrm_1 = zeros(size(Err,2),1);
            for ii = 1:size(Err,2)
                nrm_1(ii,1) = norm(Err(:,ii),inf);
            end
            fprintf('Error in J_orthogonal of U_in = %11.4e \n', max(nrm_1))
        
            Err2 = U_out' * [ U_out(dim_n+1:dim_large,:); - U_out(1:dim_n,:) ];

            nrm_2 = zeros(size(Err2,2),1);
            for ii = 1:size(Err2,2)
                nrm_2(ii,1) = norm(Err2(:,ii),inf);
            end
            fprintf('Error in J_orthogonal of U_out = %11.4e \n', max(nrm_2))
        
            Err3 = U_out' * [ U_in(dim_n+1:dim_large,:); - U_in(1:dim_n,:) ] + eye(size(U_out,2));
            nrm_3 = zeros(size(Err3,2),1);
            for ii = 1:size(Err3,2)
                nrm_3(ii,1) = norm(Err3(:,ii),inf);
            end
            fprintf('Error in J_orthogonal of U_in and U_out = %11.4e \n', max(nrm_3))
            
            Err4 = U_on_u' * [ U_on_u(dim_n+1:dim_large,:); -U_on_u(1:dim_n,:)] - ...
                [zeros(mm_unit_ew), eye(mm_unit_ew); -eye(mm_unit_ew), zeros(mm_unit_ew)];
            nrm_4 = zeros(size(Err4,2),1);
            for ii = 1:size(Err4,2)
                nrm_4(ii,1) = norm(Err4(:,ii),inf);
            end
            fprintf('Error in J_orthogonal of U_on_u = %11.4e \n', max(nrm_4))
            
            mm   = size(U_eig_space,2)/2;
            Err5 = U_eig_space' * [ U_eig_space(dim_n+1:dim_large,:); -U_eig_space(1:dim_n,:)] - ...
                [zeros(mm), eye(mm); -eye(mm), zeros(mm)];
            nrm_5 = zeros(size(Err5,2),1);
            for ii = 1:size(Err5,2)
                nrm_5(ii,1) = norm(Err5(:,ii),inf);
            end
            fprintf('Error in J_orthogonal of Eig_Space = %11.4e \n', max(nrm_5))
            
        end
        
        X_init      = randn(dim_large, 2*(dim_n - no_cong_ew));
        mtx_Y       = U_eig_space' * [ X_init(dim_n+1:dim_large,:); -X_init(1:dim_n,:) ];
        X_new       = X_init - U_eig_space * [ -mtx_Y(no_cong_ew+1:2*no_cong_ew,:); mtx_Y(1:no_cong_ew,:) ];
               
        % 
        % J_orthogonization
        %
        J_orth = 'no';
        if ( strcmp(J_orth,'yes') )
            mm = size(X_new,2) / 2;
            for ii = 1:mm
                for jj = 1:ii-1
                    
                    alpha       = X_new(:,mm+jj)'*[ X_new(dim_n+1:dim_large,ii); -X_new(1:dim_n,ii)];
                    beta        = X_new(:,jj)'*[ -X_new(dim_n+1:dim_large,ii); X_new(1:dim_n,ii)];
                    X_new(:,ii) = X_new(:,ii) + alpha * X_new(:,jj) + beta * X_new(:,mm+jj);

                    alpha          = X_new(:,mm+jj)'*[ X_new(dim_n+1:dim_large,mm+ii); -X_new(1:dim_n,mm+ii)];
                    beta           = X_new(:,jj)'*[ -X_new(dim_n+1:dim_large,mm+ii); X_new(1:dim_n,mm+ii)];
                    X_new(:,mm+ii) = X_new(:,mm+ii) + alpha * X_new(:,jj) + beta * X_new(:,mm+jj); 
                    
                end
                X_new(:,ii)    = X_new(:,ii) / norm(X_new(:,ii));
                X_new(:,mm+ii) = X_new(:,mm+ii) / (X_new(:,ii)' * [ X_new(dim_n+1:dim_large,ii+mm); -X_new(1:dim_n,ii+mm) ]);
            end
        end
        
        RR_A        = X_new' * (mtx_M * X_new);
        RR_B        = X_new' * (mtx_L * X_new);
        
        no_loss_ew      = dim_n - (mm_real + 2 * mm_cmp + mm_unit_ew);
        loss_ew         = zeros(no_loss_ew,1);
        loss_ev         = zeros(2*dim_n, no_loss_ew);
        loss_ev_out     = zeros(2*dim_n, no_loss_ew);
        ii_ew           = 0;
        
        [RR_EV0,RR_ew0] = eig(RR_A, RR_B, 'qz');
        RR_ew0          = diag(RR_ew0); 
        %
        % Find the losing eigenvalues
        %
        tol             = 1.0e-6;
        kk              = 1;
        flag_find_ew    = 1;
        while ( kk <= 10 && flag_find_ew == 1 )
            idx_ew          = find(abs(RR_ew0) <= 1+1.0e-9 & imag(RR_ew0) > -tol);
            tmp_ew          = RR_ew0(idx_ew);
            
            idx_loss_unit_ew = find(abs(abs(tmp_ew)-1) < 1.0e-10);
            idx_loss_cmp_ew  = find(abs(imag(tmp_ew)) > 1.0e-9 & abs(abs(tmp_ew)-1) >= 1.0e-10);
            idx_loss_real_ew = setdiff(1:length(tmp_ew), [idx_loss_cmp_ew; idx_loss_unit_ew]);
            mm_loss_real     = length(idx_loss_real_ew);
            mm_loss_cmp      = length(idx_loss_cmp_ew);
            mm_loss_unit_ew  = length(idx_loss_unit_ew);
            T_mm_loss        = mm_loss_real + 2 * mm_loss_cmp + mm_loss_unit_ew;
        
            if ( T_mm_loss == no_loss_ew )
                flag_find_ew = 0;
            elseif ( T_mm_loss > no_loss_ew )
                flag_find_ew = 1;
                tol          = tol / 2;
            else
                flag_find_ew = 1;
                tol          = tol * 2;
            end
            kk = kk + 1; 
        end
        RR_ew           = RR_ew0(idx_ew);
        RR_EV           = RR_EV0(:,idx_ew);
        
        %
        % Sorting eigenvalues according the angles
        %
        angle_ew = angle(RR_ew);
        idx_tmp  = find( abs(angle_ew) <= 1.0e-10 );
        if ( ~isempty(idx_tmp) )
            angle_ew(idx_tmp) = 0;
        end
        angle_ew     = mod(angle_ew,2*pi);
        [~,idx_sort] = sort(angle_ew);
        angle_ew     = angle_ew(idx_sort);
        RR_ew        = RR_ew(idx_sort);
        RR_EV        = RR_EV(:,idx_sort);
        undo_RR_ew   = zeros(length(idx_sort),1);
        no_undo_ew   = 0;
        RR_ew 
        
%         jk = 1;
%         while ( jk <= 20 && ~isempty(RR_ew) )
%             jk = jk + 1;
        kk        = 1;
        ratio     = 0.1; %0.2;
        T_mm_loss = 0;
        
        while ( kk <= 10 && T_mm_loss < no_loss_ew)
        %while ( kk <= 10 && ~isempty(RR_ew) && T_mm_loss < no_loss_ew)
                ii_ew_old    = ii_ew;
                leng_ew      = length(RR_ew);
                idx_tmp      = find(abs(angle_ew(1) - angle_ew(2:leng_ew)) < 15/180*pi);
                if ( ~isempty(idx_tmp) )
                    ratio
                    idx_tmp       = idx_tmp + 1;
                    idx_tmp2      = abs(RR_ew(1)- RR_ew(idx_tmp)) < ratio * abs(RR_ew(1));
                    idx_tmp       = idx_tmp(idx_tmp2);
                    target_ew     = [RR_ew(1); RR_ew(idx_tmp)];
                    idx_remainder = setdiff(1:leng_ew, [1; idx_tmp]);
                else
                    target_ew     = RR_ew(1);
                    idx_remainder = 2:leng_ew;
                end
         
                m_target_ew = length(target_ew);
                % Compute the target eigenpairs
                %
%                 if ( isempty(target_ew) ) 
%                     tt        = (mtx_M - target_ew(1) * mtx_L) \ (X_new * RR_EV(:,1));
%                     tt        = tt / norm(tt);
%                     Mt        = mtx_M * tt;
%                     Lt        = mtx_L * tt; 
%                     RitzValue = (tt' * Mt) / (tt' * Lt);
%                     rsdl      = norm(Mt - RitzValue * Lt);
%                     
%                     fprintf('rsdl(%24.16e+(%24.16e)1i) = %10.4e \n', real(RitzValue), ...
%                     imag(RitzValue), rsdl);
%                 
%                     loss_ew(ii_ew+1,1) = RitzValue;
%                     loss_ev(:,ii_ew+1) = tt;
%                     ii_ew              = ii_ew + 1;
%                     
%                 else 
                    if ( m_target_ew == 1 )
                        shift_S_Sinv  = sum(target_ew*0.999 + 1./(target_ew*0.999)) / length(target_ew);
                    else
                        shift_S_Sinv  = sum(target_ew + 1./target_ew) / length(target_ew);
                    end
                    tmp_sqrt_root = sqrt(shift_S_Sinv^2 - 4);
                    tmp           = [ (shift_S_Sinv+tmp_sqrt_root)/2; (shift_S_Sinv-tmp_sqrt_root)/2 ];
                    [~,idx_min]   = min(abs(tmp - RR_ew(1)));
                    new_shift     = tmp(idx_min);
                    conv_radius   = max(abs(target_ew + 1./target_ew - new_shift - 1 / new_shift));
                    idx_in_cn_reg = find( abs(ew_inside_c + 1./ew_inside_c - new_shift - 1 / new_shift) <= conv_radius );
                    eigenwanted   = length(target_ew) + min(length(idx_in_cn_reg), 30) + 6;
%                     eigenwanted   = min(2*length(target_ew),length(target_ew) + length(idx_in_cn_reg)) + 5;
                
                    W                  = mtx_NME.mtx_A - new_shift*mtx_NME.mtx_Q + (new_shift^2)*(mtx_NME.mtx_A.');
                    LU_W.Perm_amd_vec  = amd(W);
                    W_reorder          = W(LU_W.Perm_amd_vec, LU_W.Perm_amd_vec);
                    LU_W.Perm_amd      = sparse(LU_W.Perm_amd_vec,1:dim_n,ones(dim_n,1));
        
                    [LU_W.Low_L,LU_W.upper_U,LU_W.Perm_LU] = lu(W_reorder);
        
                    LU_W.Perm_LU_Tran  = LU_W.Perm_LU.';
                    LU_W.Low_L_Tran    = LU_W.Low_L.';
                    LU_W.upper_U_Tran  = LU_W.upper_U.'; 
                    LU_W.sigma         = new_shift;
        
                    if ( abs(new_shift) <= 3.0e-3 )
                        LU_W.mtx_Qp  = W;
                        LU_W.mtx_QpT = W.';
                    end
                
                    maxit_GTSHIRA = 30;
                    restart_info  = 'no';
                
                    [ ew_new, ev_new, outer_it ] = Drive_GTSHIRA( dim_n, mtx_NME.mtx_A, mtx_NME.mtx_Q, mtx_NME, ...
                        @(x,transp_flag)solve_LS_PQEP(x,transp_flag,LU_W), new_shift, eigenwanted, tolerance, maxit_GTSHIRA, ...
                        restart_info, [], [], [], [] );
                  
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
                    tol        = 1.0e-7;
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
                
                        for i1 = 1:length(idx_inside)
                            if (min(abs(lambda_in(idx_inside(i1))-ew_inside_c))/abs(lambda_in(idx_inside(i1))) < tol )
                                idx_deupli(idx_inside(i1),1) = 1;
                            end
                        end
                
                    end
                
                %end
        
                    idxx = find(idx_new_ui_ew == 1);
                    if ( ~isempty(idxx) )
                        mm3                                      = length(idxx);  
                        loss_ew(ii_ew+1:ii_ew+mm3,1)             = lambda_in(idxx);
                        loss_ev(:,ii_ew+1:ii_ew+mm3)             = ev_in_c(:,idxx); 
                        ii_ew                                    = ii_ew + mm3; 
                    
                        ew_unit_c(no_unit_ew+1:no_unit_ew+mm3,1) = lambda_in(idxx);
                        ew_unit_c(no_unit_ew+1:no_unit_ew+mm3,2) = ew_new(2*idxx);
                        ev_on_c(:,no_unit_ew+1:no_unit_ew+mm3,1) = ev_in_c(:,idxx);
                        ev_on_c(:,no_unit_ew+1:no_unit_ew+mm3,2) = ev_out_c(:,idxx);
                        no_unit_ew                               = no_unit_ew + mm3;
                
                    end
            
                    idxx              = find(idx_deupli == 0); 
                    
                    if ( isempty(idxx) )
                        fprintf('All eigenvalues are duplicated \n');
                        ratio                                             = ratio / 2;
                        undo_RR_ew(no_undo_ew+1:no_undo_ew+m_target_ew,1) = target_ew;
                        no_undo_ew                                        = no_undo_ew+m_target_ew;
                        
                        RR_ew         = RR_ew(idx_remainder);
                        if ( ~isempty(RR_ew) )
                            %RR_EV         = RR_EV(:,idx_remainder);
                            angle_ew      = angle_ew(idx_remainder);
                        else
                            RR_ew         = undo_RR_ew(1:no_undo_ew); 
                            angle_ew      = angle(RR_ew);
                            idx_tmp       = find( abs(angle_ew) <= 1.0e-10 );
                            if ( ~isempty(idx_tmp) )
                                angle_ew(idx_tmp) = 0;
                            end
                            angle_ew     = mod(angle_ew,2*pi);
                            [~,idx_sort] = sort(angle_ew);
                            angle_ew     = angle_ew(idx_sort);
                            RR_ew        = RR_ew(idx_sort);
                            ratio        = ratio / 2;
                            
                            undo_RR_ew(1:no_undo_ew) = 0;
                            no_undo_ew               = 0;
                        end
                
                    else
                
                        idx1              = find(imag(lambda_in(idxx)) > -1.0e-10);
                        
                        if ( ~isempty(idx1) ) 
                        mm2               = length(idx1);
                        idxx              = idxx(idx1);
                    %
                        idx2              = zeros(mm2*2,1); 
                        idx2(1:2:2*mm2,1) = 2*idxx-1;
                        idx2(2:2:2*mm2,1) = 2*idxx;
                        ew_new            = ew_new(idx2);
                        ev_new            = ev_new(:,idx2);
            
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
                    
                            loss_ew(ii_ew+1:ii_ew+mm2,1)     = lambda_in(idxx);
                            loss_ev(:,ii_ew+1:ii_ew+mm2)     = ev_in_c(:,idxx);
                            loss_ev_out(:,ii_ew+1:ii_ew+mm2) = ev_out_c(:,idxx);
                            ii_ew                            = ii_ew + mm2; 
                    
                            ew_inside_c(no_ew_in_c+1:no_ew_in_c+mm2,1) = lambda_in(idxx);
                            %ev_inside_c(:,no_ew_in_c+1:no_ew_in_c+mm2) = ev_in_c(:,idxx);
                            no_ew_in_c                                 = no_ew_in_c + mm2;  
                    
                        elseif ( leng_idx_rsdl > 0 )
                            loss_ew(ii_ew+1:ii_ew+leng_idx_rsdl,1)     = lambda_in(idxx(idx_rsdl));
                            loss_ev(:,ii_ew+1:ii_ew+leng_idx_rsdl)     = ev_in_c(:,idxx(idx_rsdl));
                            loss_ev_out(:,ii_ew+1:ii_ew+leng_idx_rsdl) = ev_out_c(:,idxx(idx_rsdl));
                            ii_ew                                      = ii_ew + leng_idx_rsdl; 
                    
                            ew_inside_c(no_ew_in_c+1:no_ew_in_c+leng_idx_rsdl,1) = lambda_in(idxx(idx_rsdl)); 
                            no_ew_in_c                                           = no_ew_in_c + leng_idx_rsdl; 
                    
                        end
                                
                        no_computed_ew = ii_ew - ii_ew_old;
                        if ( no_computed_ew ~= length(target_ew) )
%                     if ( no_computed_ew > length(target_ew) )
%                         fprintf('Error for finding losing eigenvalues \n');
%                         [leng_idx_rsdl ii_ew ii_ew_old no_computed_ew length(target_ew)]
%                         X_sol = [];
%                         return
%                     else
                            idx_computed_ew = zeros(length(RR_ew),1);
                            tmp             = loss_ew(ii_ew_old+1:ii_ew_old+no_computed_ew);
                            avg_tmp         = sum(tmp) / no_computed_ew;
                            [~,idx_s]       = sort(abs(tmp - avg_tmp));
                            tmp             = tmp(idx_s);
                            for jj = 1:no_computed_ew
                                [min_val,idx_min] = min(abs(RR_ew - tmp(jj)));
                                if ( min_val / abs(tmp(jj)) <= 1.0e-2 )
%                             [min_val,idx_min] = min(abs(RR_ew - loss_ew(ii_ew_old+jj)));
%                             if ( min_val / abs(loss_ew(ii_ew_old+jj)) <= 1.0e-2 )
                                    idx_computed_ew(idx_min,1) = 1;
                                else
                                    fprintf('Error in form_eigspace \n');
                                    [loss_ew(ii_ew_old+jj) RR_ew(idx_min,1)  min_val/abs(tmp(jj))]
                                end
                            end
                            idx_remainder = find( idx_computed_ew == 0 );
%                     end
                        end
                        end
                        
                        RR_ew         = RR_ew(idx_remainder);
                        %RR_EV         = RR_EV(:,idx_remainder);
                        angle_ew      = angle_ew(idx_remainder);
                
                        idx_finding_unit_ew = find(abs(abs(loss_ew(1:ii_ew))-1) < 1.0e-10);
                        idx_finding_cmp_ew  = find(abs(imag(loss_ew(1:ii_ew))) > 1.0e-9 & abs(abs(loss_ew(1:ii_ew))-1) >= 1.0e-10);
                        idx_finding_real_ew = setdiff(1:ii_ew, [idx_finding_cmp_ew; idx_finding_unit_ew]);
                        mm_finding_real     = length(idx_finding_real_ew);
                        mm_finding_cmp      = length(idx_finding_cmp_ew);
                        mm_finding_unit_ew  = length(idx_finding_unit_ew);
                        T_mm_loss           = mm_finding_real + 2 * mm_finding_cmp + mm_finding_unit_ew;
                        
                        if ( isempty(RR_ew) && no_undo_ew > 0 && T_mm_loss < no_loss_ew) 
                            RR_ew         = undo_RR_ew(1:no_undo_ew); 
                            angle_ew      = angle(RR_ew);
                            idx_tmp       = find( abs(angle_ew) <= 1.0e-10 );
                            if ( ~isempty(idx_tmp) )
                                angle_ew(idx_tmp) = 0;
                            end
                            angle_ew     = mod(angle_ew,2*pi);
                            [~,idx_sort] = sort(angle_ew);
                            angle_ew     = angle_ew(idx_sort);
                            RR_ew        = RR_ew(idx_sort);
                            ratio        = ratio / 2;
                            
                            undo_RR_ew(1:no_undo_ew) = 0;
                            no_undo_ew               = 0;
                        end
                
                    end
                %end
                
            %idx_remainder = setdiff(1:length(RR_ew), [1; idx_tmp+1]);
                
                kk            = kk + 1;
                
        end
        %end
        
%         idx_finding_unit_ew = find(abs(abs(loss_ew(1:ii_ew))-1) < 1.0e-10);
%         idx_finding_cmp_ew  = find(abs(imag(loss_ew(1:ii_ew))) > 1.0e-9 & abs(abs(loss_ew(1:ii_ew))-1) >= 1.0e-10);
%         idx_finding_real_ew = setdiff(1:ii_ew, [idx_finding_cmp_ew; idx_finding_unit_ew]);
%         mm_finding_real     = length(idx_finding_real_ew);
%         mm_finding_cmp      = length(idx_finding_cmp_ew);
%         mm_finding_unit_ew  = length(idx_finding_unit_ew);
        U_in_finding        = zeros(dim_large,mm_finding_real+2*mm_finding_cmp);
        U_out_finding       = zeros(dim_large,2);
        
        U_in_finding (:,1:mm_finding_real) = real(loss_ev(:,idx_finding_real_ew)); 
        for ii = 1:mm_finding_real
            U_in_finding(:,ii)  = U_in_finding(:,ii) / norm(U_in_finding(:,ii)); 
        end
        kk                         = mm_finding_real;
        
        for ii = 1:mm_finding_cmp
            [ U_in_finding(:,kk+1:kk+2), U_out_finding(:,kk+1:kk+2) ] = J_orthogonal( [real(loss_ev(:,idx_finding_cmp_ew(ii))) ...
                imag(loss_ev(:,idx_finding_cmp_ew(ii)))], [real(loss_ev_out(:,idx_finding_cmp_ew(ii))) ...
                imag(loss_ev_out(:,idx_finding_cmp_ew(ii)))], dim_n );
            kk                                                        = kk + 2;
        end
        
        no_cong_ew = no_cong_ew + mm_finding_real + 2 * mm_finding_cmp + mm_finding_unit_ew;
        if ( no_cong_ew == dim_n )
            find_all_eigenvalue = 'yes';
            U_in                = [ U_in, U_in_finding];
        end
        
    end
    
    if ( strcmp(find_all_eigenvalue,'yes' ) )
        
        tmp_vec2 = zeros(size(ew_unit_c,1),2);
        for j2 = 1:size(ew_unit_c,1)
            for i2 = 1:2
                tmp_vec2(j2,i2) = 1i * (ev_on_c(1:dim_n,j2,i2)'*(2 * ew_unit_c(j2,i2) * ...
                     (mtx_NME.mtx_A.'*ev_on_c(1:dim_n,j2,i2)) - mtx_NME.mtx_Q*ev_on_c(1:dim_n,j2,i2))); 
            end
        end
        
        idx = find( abs(imag(tmp_vec2)) > 1.0e-7 );
        if ( ~isempty(idx) )
            fprintf('Error in choosing unimodular eigenvalues with index \n');
            display(idx)
        else
            idx1           = find(real(tmp_vec2(:,1)) > 0 & real(tmp_vec2(:,2)) < 0);
            idx2           = find(real(tmp_vec2(:,2)) > 0 & real(tmp_vec2(:,1)) < 0);
            if ( length(idx1) + length(idx2) ~= size(ew_unit_c,1))
                fprintf('Error in choosing unimodular eigenvalues \n');
            else
                target_ev_on_c = zeros(2*dim_n,size(ew_unit_c,1));
                target_unit_ew = [ ew_unit_c(idx1,1); ew_unit_c(idx2,2) ];
                target_ev_on_c(:,1             :length(idx1)             ) = ev_on_c(:,idx1,1);
                target_ev_on_c(:,1+length(idx1):length(idx1)+length(idx2)) = ev_on_c(:,idx2,2); 
            end
        end

        X_subspace = [ U_in target_ev_on_c ]; 
        X_sol      = X_subspace(dim_n+1:dim_large,:) / X_subspace(1:dim_n,:);
        X_sol      = ( X_sol + X_sol.') / 2;
        Err        = mtx_NME.mtx_A.' * (X_sol \ mtx_NME.mtx_A);
        nrm_2      = norm(Err,inf);
        Err        = Err + X_sol - mtx_NME.mtx_Q;
        RelErr     = norm(Err,inf) / (norm(X_sol,inf) + nrm_2 + norm(mtx_NME.mtx_Q,inf));
    else
        fprintf('Not all eigenvectors are computed \n');
        loss_ew
        RR_ew
        X_sol  = [];
        RelErr = 1;
    end
            
end

