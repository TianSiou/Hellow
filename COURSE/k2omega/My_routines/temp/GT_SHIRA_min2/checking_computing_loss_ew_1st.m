function [ ew_new_1, ev_new, rsdl, cong_flag, ew_unit_c, ev_on_c, ew_inside_c_new, ev_inside_c_new, ...
         no_ew_in_c, next_ref_ev_outside_c, no_ref_ev, flag_conv, arrive_left_tag, no_unit_ew ] = checking_computing_loss_ew_1st( dim_n, ...
         mtx_NME, lambda_in, prev_ew, ew_new, ev_new, ref_ew_inside_c, sigma, eigenwanted, ...
         tolerance, conv_radius, rad, rad_old, ew_unit_c, no_ew_in_c, ew_inside_c_new, ...
         ev_on_c, ev_inside_c_new, next_ref_ev_outside_c, no_ref_ev, mtx_M, mtx_L, ew_no_dupli, ...
         deflation_enforce, min_rad_neg_real_ew, arrive_left_tag, no_unit_ew )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

flag_conv   = 0;
cong_flag   = 0;

idx1        = find( abs(lambda_in) < rad_old, 1 );
if ( isempty(idx1) )
    ew_new_1 = [];
    rsdl     = 0;
    ev_new   = [];
    return
end

tol_angle   = 1.0e-10;

idx_angle        = find(imag(lambda_in) >= -1.0e-10 & abs(lambda_in) <= 1);
left_ang_thm_pt  = esstimating_fm_in_conv_circle(conv_radius, sigma, lambda_in(idx_angle), 'left');
right_ang_thm_pt = esstimating_fm_in_conv_circle(conv_radius, sigma, lambda_in(idx_angle), 'right');

t1           = angle(sigma);
if ( abs(t1) <= tol_angle )
    t1 = 0;
end
    
if ( ~isempty(idx_angle) )
    tmp_ang2              = angle(lambda_in(idx_angle));
    idxx                  =  abs(tmp_ang2) <= tol_angle ;
    tmp_ang2(idxx)        = 0;
    tmp_ang2              = mod(tmp_ang2,2*pi); 
    left_angle            = max(0.95 * left_ang_thm_pt, max(tmp_ang2)); %0.95 * max(left_ang_thm_pt, max(tmp_ang2));
    [min_ang,i_right_ang] = min(tmp_ang2);
    right_angle           = min(right_ang_thm_pt, min_ang); 
else
    left_angle  = 0.95 * left_ang_thm_pt;
    right_angle = right_ang_thm_pt;
end

ref_angle = angle(ref_ew_inside_c);
idxx      = find( abs(ref_angle) <= tol_angle );
if ( ~isempty(idxx) )
    ref_angle(idxx) = 0;
end
ref_angle  = mod(ref_angle,2*pi);

%
%
%
[ doing_GTSHIRA_again, new_shift ] = checking_loss_ew_1st( lambda_in, prev_ew, ref_ew_inside_c, ...
    sigma, conv_radius, rad, rad_old, ew_no_dupli );

if ( strcmp(doing_GTSHIRA_again,'yes') )
    fprintf('\n Need to re-do GTSHIRA with shift = %11.4e + (%11.4e) * 1i \n\n',real(new_shift), imag(new_shift));
    
     maxit              = 70;

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
            
    nmax               = min(150, maxit + round(0.9 / abs(new_shift)));

    [ ew_T, residual, R, H, Y, Z, nrm_Err_N ] = estimate_by_GTSHIRA( mtx_NME.mtx_A, mtx_NME.mtx_Q, ...
                new_shift, mtx_NME, @(x,transp_flag)solve_LS_PQEP(x,transp_flag,LU_W), dim_n, nmax );
    
    idx_pre_ew      = find( abs(pre_ew) <= 0.99 );
    min_right_real  = min(real(pre_ew(idx_pre_ew)));
    idx_ew          =  abs(lambda_in) <= 0.99 ;
    max_left_real   = max(real(lambda_in(idx_ew)));
    
    mu0             = new_shift + 1 / new_shift;
    dist_right_ew   = min(abs(prev_ew + 1./prev_ew - mu0));
    dist_left_ew    = min(abs(lambda_in + 1./lambda_in - mu0));
    est_conv_radius = max(abs(ew_T + 1./ew_T - mu0));
     
    if ( est_conv_radius >= max(dist_right_ew,dist_left_ew)*1.02 )
        idx_target_ew = find( real(ew_T) >= max_left_real-0.05 & real(ew_T) <= min_right_real+0.05);
        eigenwanted   = length(idx_target_ew) + 3;
    else
        eigenwanted   = length(ew_T) + 15;
    end
    
    restart_info    = 'yes';
    maxit_GTSHIRA   = 30;
       
    fprintf('GTSHIRA with sigma = %11.4e + (%11.4e)*1i \n', real(new_shift), imag(new_shift));
    [ ew_new, ev_new, outer_it ] = Drive_GTSHIRA( dim_n, mtx_NME.mtx_A, mtx_NME.mtx_Q, mtx_NME, ...
         @(x,transp_flag)solve_LS_PQEP(x,transp_flag,LU_W), new_shift, eigenwanted, tolerance, ...
                        maxit_GTSHIRA, restart_info, R, H, Y, Z );
    
    if ( ~isempty(ew_new) ) % GTSHIRA converges or not
        cong_flag     = 0; 
                    
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
        end
            
        idxx              = find(idx_deupli == 0); 
        mm2               = length(idxx);
        percentage        = (leng_ew/2-mm2)/(leng_ew/2); 
        fprintf('Re_GTSHIRA, percentage of deuplicated eigenvalues = %11.4e \n', percentage)
 
        idx1              = find(imag(lambda_in(idxx)) > -1.0e-10);
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
    % ==================================================                

        
end

