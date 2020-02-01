function [ ew_new_1, ev_new, rsdl, cong_flag, ew_unit_c, ev_on_c, ew_inside_c_new, ev_inside_c_new, ...
         no_ew_in_c, next_ref_ev_outside_c, no_ref_ev, flag_conv, arrive_left_tag, no_unit_ew ] = checking_computing_loss_ew( dim_n, ...
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

% mu_0        = sigma + 1 / sigma;
% 
% tmp         = [ (-(conv_radius-mu_0)+sqrt((conv_radius-mu_0)^2-4))/(2*rad); 
%                 (-(conv_radius-mu_0)-sqrt((conv_radius-mu_0)^2-4))/(2*rad);
%                 ((conv_radius+mu_0)+sqrt((conv_radius+mu_0)^2-4))/(2*rad);
%                 ((conv_radius+mu_0)-sqrt((conv_radius+mu_0)^2-4))/(2*rad) ];
% idx_tmp     = imag(tmp) >= -1.0e-10;
% tmp_ang1    = angle(tmp(idx_tmp));
% idxx        = find( abs(tmp_ang1) <= tol_angle );
% if ( ~isempty(idxx) )
%     tmp_ang1(idxx) = 0;
% end
% tmp_ang1    = mod(tmp_ang1,2*pi);

%idx_angle   = find(imag(lambda_in) >= -1.0e-10 & abs(lambda_in) <= 1); 

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

idx2       =  abs(ref_ew_inside_c) <= rad_old & ref_angle >= right_angle & ref_angle <= left_angle;

top_ref_ew = ref_ew_inside_c(idx2);
%
%
%
[ doing_GTSHIRA_again, new_shift ] = checking_loss_ew( lambda_in, prev_ew, ref_ew_inside_c, ...
    sigma, conv_radius, rad, rad_old, ew_no_dupli );

rsdl_tol = 1.0e-9;

if ( strcmp(doing_GTSHIRA_again,'yes') )
    fprintf('\n Need to re-do GTSHIRA with shift = %11.4e + (%11.4e) * 1i \n\n',real(new_shift), imag(new_shift));
                    
    eigenwanted = max(15,eigenwanted);
    
    [ deflation_enforce_1, lambda_in1, ew_new_1, ev_new_1, conv_radius_1, rsdl, no_deuplicate_ew_1, ...
        cong_flag, ew_unit_c, ev_on_c, ew_inside_c_new, ev_inside_c_new, ...
        no_ew_in_c, next_ref_ev_outside_c, no_ref_ev, idxx, no_unit_ew ] = compute_ew_near_sigma( dim_n, mtx_NME, ...
        new_shift, deflation_enforce, rad, rad_old, ew_new, ev_new, eigenwanted, tolerance, ew_unit_c, ...
        no_ew_in_c, ew_inside_c_new, ev_on_c, ev_inside_c_new, next_ref_ev_outside_c, ...
        no_ref_ev, mtx_M, mtx_L, no_unit_ew );
                    
    if ( max(rsdl) > rsdl_tol )
        flag_conv = 1;
        return
    end
    
    if ( ~isempty(ew_new_1) )
                
        angle_new_ew = angle(lambda_in1);
        idx_ang      = abs(angle_new_ew) <= tol_angle;
        if ( ~isempty(idx_ang) )
            angle_new_ew(idx_ang) = 0;
        end
        angle_new_ew = mod(angle_new_ew, 2*pi);
    
        %
        % Arrive eigenvalues in rad_old or not
        %
        idx_top = find(abs(top_ref_ew + 1./top_ref_ew - new_shift - 1/new_shift) <= conv_radius_1, 1); 
        if ( isempty(idx_top) && ~isempty(top_ref_ew) )
            idx_new_ew  = find(abs(ew_new_1) <= 1);
            if ( ~isempty(idx_new_ew) )
                [val,ii]    = min(abs(top_ref_ew));
                shift       = (max(abs(ew_new_1(idx_new_ew))) + val) / 2;
                t1          = angle(new_shift);
                t2          = angle(top_ref_ew(ii));
                if ( abs(t1) <= tol_angle )
                    t1 = 0;
                end
                if ( abs(t2) <= tol_angle )
                    t2 = 0;
                end
                angle_shift = (min(angle_new_ew) + min(mod(t1,2*pi), mod(t2, 2*pi))) / 2; 
                shift       = shift * exp(1i * angle_shift);
        
                fprintf('\n Compute the eigenvalues between current step and previous radius with shift = %11.4e+(%11.4e)*1i \n', real(shift), imag(shift));
                
                [ deflation_enforce_1, lambda_in2, ew_new_2, ev_new_1, conv_radius_1, rsdl, no_deuplicate_ew_1, ...
                    cong_flag, ew_unit_c, ev_on_c, ew_inside_c_new, ev_inside_c_new, ...
                    no_ew_in_c, next_ref_ev_outside_c, no_ref_ev, idxx, no_unit_ew ] = compute_ew_near_sigma( dim_n, mtx_NME, ...
                    shift, 'no', rad, rad_old, ew_new, ev_new, eigenwanted, tolerance, ew_unit_c, ...
                    no_ew_in_c, ew_inside_c_new, ev_on_c, ev_inside_c_new, next_ref_ev_outside_c, ...
                    no_ref_ev, mtx_M, mtx_L, no_unit_ew );
        
                if ( max(rsdl) > rsdl_tol )
                    flag_conv = 1;
                    return
                end
        
                ew_new_1   = [ ew_new_1; ew_new_2 ];
                lambda_in1 = [ lambda_in1; lambda_in2 ];
            end
        end
     
        %
        % Arrive eigenvalues in rad or not
        %
        idx_top = find(abs(lambda_in + 1./lambda_in - new_shift - 1/new_shift) <= conv_radius_1, 1);
        if ( isempty(idx_top) ) 
            [val,ii]    = min(abs(ew_new_1));
            shift       = (max(abs(lambda_in)) + val) / 2;
            t1          = angle(new_shift);
            t2          = angle(ew_new_1(ii));
            if ( abs(t1) <= tol_angle )
                t1 = 0;
            end
            if ( abs(t2) <= tol_angle )
                t2 = 0;
            end
            angle_shift = (right_angle + min(mod(t1,2*pi), mod(t2, 2*pi))) / 2;
            shift       = shift * exp(1i * angle_shift);
        
            fprintf('\n Compute the eigenvalues between current step and radius with shift = %11.4e+(%11.4e)*1i \n', real(shift), imag(shift));
            
           [ deflation_enforce_1, lambda_in2, ew_new_2, ev_new_1, conv_radius_1, rsdl, no_deuplicate_ew_1, ...
                cong_flag, ew_unit_c, ev_on_c, ew_inside_c_new, ev_inside_c_new, ...
                no_ew_in_c, next_ref_ev_outside_c, no_ref_ev, idxx, no_unit_ew ] = compute_ew_near_sigma( dim_n, mtx_NME, ...
                shift, 'no', rad, rad_old, ew_new, ev_new, eigenwanted, tolerance, ew_unit_c, ...
                no_ew_in_c, ew_inside_c_new, ev_on_c, ev_inside_c_new, next_ref_ev_outside_c, ...
                no_ref_ev, mtx_M, mtx_L, no_unit_ew );
        
            if ( max(rsdl) > rsdl_tol )
                flag_conv = 1;
                return
            end
        
            ew_new_1   = [ ew_new_1; ew_new_2 ];
            lambda_in1 = [ lambda_in1; lambda_in2 ];
        end
        
    else
        lambda_in1 = [];
    end
else
    
%     t1           = angle(sigma);
%     if ( abs(t1) <= tol_angle )
%         t1 = 0;
%     end
    %
    % Check the eigenvalues in right top corner are convergent or not
    %
    % Redefine the region of the angle from right_angle_new to left_ang_new
    %
    left_ang_new = (right_angle + mod(t1,2*pi)) / 2;
    
    if ( ~isempty(ew_no_dupli) )
            %idx_tmp1     =  imag(prev_ew) >= - tol_angle ;
            angle_pre_ew = angle(ew_no_dupli); %angle(prev_ew(idx_tmp1));
            idx_tmp1     = find( abs(angle_pre_ew) <= tol_angle );
            if ( ~isempty(angle_pre_ew(idx_tmp1)) )
                angle_pre_ew(idx_tmp1) = 0;
            end
            angle_pre_ew    = mod(angle_pre_ew, 2*pi);
            right_angle_new = min(angle_pre_ew); 
            %angle_pre_ew    = mod(angle_pre_ew, 2*pi);
            %right_angle_new = min(angle_pre_ew); %max(right_angle, max(angle_pre_ew))
    else
            right_angle_new = right_angle;
    end
        
    idx2         = abs(ref_ew_inside_c) <= rad_old & ref_angle >= right_angle_new & ref_angle <= left_ang_new;

    top_ref_ew   = ref_ew_inside_c(idx2);
    idx_top      = find(abs(top_ref_ew + 1./top_ref_ew - sigma - 1/sigma) <= conv_radius, 1);
    
    if ( isempty(idx_top) && ~isempty(top_ref_ew) )
        
        idx2        =  find(angle_pre_ew >= right_angle_new & angle_pre_ew <= left_ang_new);
        if ( ~isempty(idx2) )
            shift       = (11*min(abs(top_ref_ew)) + 9*max(abs(ew_no_dupli(idx2)))) / 20;
            angle_shift = left_ang_new; %(right_angle_new + 19*left_ang_new) / 20;
            shift       = shift * exp(1i * angle_shift); 
            
            fprintf('\n Compute the eigenvalues in right top corner with shift = %11.4e+(%11.4e)*1i \n', real(shift), imag(shift));
            
            [ deflation_enforce_1, lambda_in1, ew_new_1, ev_new1, conv_radius_1, rsdl, no_deuplicate_ew_1, ...
                cong_flag, ew_unit_c, ev_on_c, ew_inside_c_new, ev_inside_c_new, ...
                no_ew_in_c, next_ref_ev_outside_c, no_ref_ev, idxx, no_unit_ew ] = compute_ew_near_sigma( dim_n, mtx_NME, ...
                shift, 'no', rad, rad_old, ew_new, ev_new, 10, tolerance, ew_unit_c, ...
                no_ew_in_c, ew_inside_c_new, ev_on_c, ev_inside_c_new, next_ref_ev_outside_c, ...
                no_ref_ev, mtx_M, mtx_L, no_unit_ew );
        
            if ( max(rsdl) > rsdl_tol )
                flag_conv = 1;
                return
            end
        else
            ew_new_1   = [];
            rsdl       = 0;
            ev_new     = [];
            lambda_in1 = [];
        end
        
    else
        
        ew_new_1   = [];
        rsdl       = 0;
        ev_new     = [];
        lambda_in1 = [];
    end
    
    %
    % Check the eigenvalues in left top corner are convergent or not
    %
    % Redefine the region of the angle from right_angle_new to left_ang_new
    %
    
    left_angle    = left_angle / 0.95;
    right_ang_new = (left_angle + mod(t1,2*pi)) / 2;
    idx2          =  abs(ref_ew_inside_c) <= rad_old & ref_angle >= right_ang_new & ref_angle <= left_angle;
    top_ref_ew    = ref_ew_inside_c(idx2);
    
    if ( ~isempty(idx2) )
        idx_top   = find(abs(top_ref_ew + 1./top_ref_ew - sigma - 1/sigma) <= conv_radius, 1);
    else
        idx_top   = [];
    end

    if ( isempty(idx_top) && ~isempty(top_ref_ew) )
        
        idx2        =  find(tmp_ang2 >= right_ang_new & tmp_ang2 <= left_angle);
        if ( ~isempty(idx2) ) 
            shift       = (min(abs(top_ref_ew)) + max(abs(lambda_in(idx_angle(idx2))))) / 2; 
            angle_shift = min(178.0/180*pi, (2*right_ang_new + left_angle) / 3);
            shift       = shift * exp(1i * angle_shift); 
        
            fprintf('\n Compute the eigenvalues in left top corner with shift = %11.4e+(%11.4e)*1i \n', real(shift), imag(shift));
            
            [ deflation_enforce_1, lambda_in2, ew_new_2, ev_new1, conv_radius_1, rsdl, no_deuplicate_ew_1, ...
                cong_flag, ew_unit_c, ev_on_c, ew_inside_c_new, ev_inside_c_new, ...
                no_ew_in_c, next_ref_ev_outside_c, no_ref_ev, idxx, no_unit_ew ] = compute_ew_near_sigma( dim_n, mtx_NME, ...
                shift, 'no', rad, rad_old, ew_new, ev_new, 5, tolerance, ew_unit_c, ...
                no_ew_in_c, ew_inside_c_new, ev_on_c, ev_inside_c_new, next_ref_ev_outside_c, ...
                no_ref_ev, mtx_M, mtx_L, no_unit_ew );
        
            if ( ~isempty(ew_new_2) )
                ew_new_1   = [ ew_new_1; ew_new_2 ];
                lambda_in1 = [ lambda_in1; lambda_in2 ];
            
                if ( max(rsdl) > rsdl_tol )
                    flag_conv = 1;
                    return
                end
            end
        end
        
    end
    
    if ( ~isempty(find(abs(lambda_in(idx_angle)) >= rad, 1)) )
    %
    % Check the eigenvalues in right bottom corner are convergent or not
    %
    % Redefine the region of the angle from right_angle_new to left_ang_new
    %
    idx1 = find( abs(prev_ew) <= rad_old );
    if ( ~isempty(prev_ew) && ~isempty(idx1) )  
        idx3 = find(abs(prev_ew(idx1) + 1./prev_ew(idx1) - (sigma + 1 / sigma)) <= conv_radius, 1);
        if ( isempty(idx3) )
            angle_prev_ew = angle(prev_ew(idx1));
            idx_angle_prev_ew = find(abs(angle_prev_ew) <= tol_angle);
            if ( ~isempty(idx_angle_prev_ew) )
                angle_prev_ew(idx_angle_prev_ew) = 0;
            end
            angle_prev_ew = mod(angle_prev_ew,2*pi);
            [max_ang, ij] = max(angle_prev_ew); 
            shift         = (abs(prev_ew(idx1(ij))) + abs(lambda_in(idx_angle(i_right_ang)))) / 2;
            shift         = shift * exp(1i*(2*right_angle + max_ang) / 3);
            
            fprintf('\n Compute the eigenvalues in right bottomp corner with shift = %11.4e+(%11.4e)*1i \n', real(shift), imag(shift));
            
            [ deflation_enforce_1, lambda_in2, ew_new_2, ev_new1, conv_radius_1, rsdl, no_deuplicate_ew_1, ...
                cong_flag, ew_unit_c, ev_on_c, ew_inside_c_new, ev_inside_c_new, ...
                no_ew_in_c, next_ref_ev_outside_c, no_ref_ev, idxx, no_unit_ew ] = compute_ew_near_sigma( dim_n, mtx_NME, ...
                shift, 'no', rad, rad_old, ew_new, ev_new, 5, tolerance, ew_unit_c, ...
                no_ew_in_c, ew_inside_c_new, ev_on_c, ev_inside_c_new, next_ref_ev_outside_c, ...
                no_ref_ev, mtx_M, mtx_L, no_unit_ew );
        
            ew_new_1   = [ ew_new_1; ew_new_2 ];
            lambda_in1 = [ lambda_in1; lambda_in2 ];
            
            if ( max(rsdl) > rsdl_tol ) %( isempty(rsdl) || max(rsdl) > 1.0e-7 )
                flag_conv = 1;
                return
            end
        end
    end
    end
    
end
 
if ( ~isempty(lambda_in1) ) 
    if ( ~isempty(idx_angle) ) 
        idx_angle_neg_ew = find( abs(tmp_ang2 - pi) <= 0.5/180*pi );
        %
        % All eigenvalues are negative real.
        %
        if ( ~isempty(idx_angle_neg_ew) ) 
            negative_real_ew = lambda_in(idx_angle(idx_angle_neg_ew));
            angle_new_ew     = mod(angle(lambda_in1), 2*pi);
            idx_new_neg_ew   = find( abs(angle_new_ew - pi) <= 0.5/180*pi );
            if ( ~isempty(idx_new_neg_ew) ) 
                new_negative_real_ew = lambda_in1(idx_new_neg_ew);
                
                min_negative_real_ew     = min(real(negative_real_ew));
                max_negative_real_ew     = max(real(negative_real_ew));
                min_new_negative_real_ew = min(real(new_negative_real_ew));
                max_new_negative_real_ew = max(real(new_negative_real_ew));
                
                if ( max_new_negative_real_ew ~= min_new_negative_real_ew)
                    per_leng_ew = length(idx_new_neg_ew) / abs(max_new_negative_real_ew - min_new_negative_real_ew);
                else
                    per_leng_ew = 0;
                end
                
                if ( min_new_negative_real_ew <= min_negative_real_ew )
                    if ( max_new_negative_real_ew < min_negative_real_ew )
                        % Loss some eigenvalues 
                        d_ew_sigma      = abs(max_new_negative_real_ew + rad);
                        d_left_right_ew = abs(max_new_negative_real_ew - min_negative_real_ew);
%                         [max_new_negative_real_ew rad min_negative_real_ew]
%                         d_left_right_ew / d_ew_sigma
                        if ( d_left_right_ew >= 0.9 * d_ew_sigma) 
                            idx_neg_ew = find(abs(imag(lambda_in)) < 1.0e-7 & real(lambda_in) <= -rad & real(lambda_in) > -rad_old );
                            if ( isempty(idx_neg_ew) )
                                new_shift    = (0.7*max_new_negative_real_ew + 0.3*min_negative_real_ew);
                            else
                                tt           = 1 / (length(idx_neg_ew) / (d_ew_sigma / (rad_old - rad)))^(0.3); 
                                rat          = 0.5 * min(1.55, 2-tt);
                                new_shift    = (rat*max_new_negative_real_ew + (1-rat)*min_negative_real_ew);
                                fprintf('tt = %11.4e, rat = %11.4e \n', tt, rad);
                            end
                        else
                            new_shift    = (max_new_negative_real_ew + min_negative_real_ew) / 2;
                        end
                        
                        new_shift    = abs(new_shift) * exp(1i * 179.5 / 180 * pi);
                        no_wanted_ew = round(per_leng_ew * abs(max_new_negative_real_ew - min_negative_real_ew));
                        no_wanted_ew = max(no_wanted_ew, 5);
                        fprintf('\n case (1) with no_wanted_ew = %3.0f in check_computing_loss_ew \n',no_wanted_ew);
                    else
                        new_shift    = 0;
                        no_wanted_ew = 10;
                    end

                    % Remove following checking at 2016/12/3
                    %
                    if ( min_new_negative_real_ew <= min_rad_neg_real_ew && mod(angle(sigma), 2*pi) > 177/180*pi)
                        arrive_left_tag = 'yes';
                    end
                else
                    if ( max_negative_real_ew < min_new_negative_real_ew )
                        % Loss some eigenvalues 
                        new_shift    = (max_negative_real_ew + min_new_negative_real_ew) / 2;
                        new_shift    = abs(new_shift) * exp(1i * 179.5 / 180 * pi);
                        no_wanted_ew = round(per_leng_ew * abs(max_negative_real_ew - min_new_negative_real_ew));
                        no_wanted_ew = max(no_wanted_ew, 5);
                        fprintf('case (2) with no_wanted_ew = %3.0f in check_computing_loss_ew \n',no_wanted_ew);
                    else
                        new_shift    = 0;
                        no_wanted_ew = 10;
                    end
                    
%                     if ( min_negative_real_ew <= min_rad_neg_real_ew )
%                         arrive_left_tag = 'yes';
%                     end
                end
%            end
            
                if ( abs(new_shift) > 0 )
                    fprintf('\n Compute the eigenvalues in x-axis with shift = %11.4e+(%11.4e)*1i \n', real(new_shift), imag(new_shift));
                
                    [ deflation_enforce_1, lambda_in2, ew_new_2, ev_new1, conv_radius_1, rsdl, no_deuplicate_ew_1, ...
                        cong_flag, ew_unit_c, ev_on_c, ew_inside_c_new, ev_inside_c_new, ...
                        no_ew_in_c, next_ref_ev_outside_c, no_ref_ev, idxx, no_unit_ew ] = compute_ew_near_sigma( dim_n, mtx_NME, ...
                        new_shift, 'no', rad, rad_old, ew_new, ev_new, no_wanted_ew, tolerance, ew_unit_c, ...
                        no_ew_in_c, ew_inside_c_new, ev_on_c, ev_inside_c_new, next_ref_ev_outside_c, ...
                        no_ref_ev, mtx_M, mtx_L, no_unit_ew );
        
                    ew_new_1   = [ ew_new_1; ew_new_2 ]; 
            
                    if ( isempty(rsdl) || max(rsdl) > rsdl_tol )
                        flag_conv = 1; 
                    end
            
                end
            end
        end
    end
else
    if ( t1 >= 178/180*pi )
        if ( ~isempty(tmp_ang2) )
            idx_angle_neg_ew = find( abs(tmp_ang2 - pi) <= 0.5/180*pi );
            if ( ~isempty(idx_angle_neg_ew) ) 
                idx_real_ew_prev = find(abs(imag(prev_ew)) < 1.0e-7 & real(prev_ew) < 0);
                
                if ( ~isempty(idx_real_ew_prev) )
                    prev_real_ew     = prev_ew(idx_real_ew_prev);
                    if ( abs(prev_real_ew + 1./prev_real_ew - sigma - 1/sigma) > conv_radius )
                        real_neg_ew = lambda_in(idx_angle(idx_angle_neg_ew));
                        max_val_cur = max(real(real_neg_ew));
                        min_val_pre = min(real(prev_real_ew));
                        
                        if (max_val_cur < min_val_pre)
                            new_shift = abs(max_val_cur + min_val_pre) / 2 * exp(1i*t1);
                        else
                            min_val_cur = min(real(real_neg_ew));
                            max_val_pre = max(real(prev_real_ew));
                            new_shift   = abs(min_val_cur + max_val_pre) / 2 * exp(1i*t1);
                        end
                        
                        no_wanted_ew = 50;
                        
                        fprintf('\n Compute the eigenvalues in x-axis with shift = %11.4e+(%11.4e)*1i \n', real(new_shift), imag(new_shift));
                
                       [ deflation_enforce_1, lambda_in1, ew_new_1, ev_new, conv_radius_1, rsdl, no_deuplicate_ew_1, ...
                           cong_flag, ew_unit_c, ev_on_c, ew_inside_c_new, ev_inside_c_new, ...
                           no_ew_in_c, next_ref_ev_outside_c, no_ref_ev, idxx, no_unit_ew ] = compute_ew_near_sigma( dim_n, mtx_NME, ...
                           new_shift, 'no', rad, rad_old, ew_new, ev_new, no_wanted_ew, tolerance, ew_unit_c, ...
                           no_ew_in_c, ew_inside_c_new, ev_on_c, ev_inside_c_new, next_ref_ev_outside_c, ...
                           no_ref_ev, mtx_M, mtx_L, no_unit_ew );
                    
                       if ( isempty(rsdl) || max(rsdl) > rsdl_tol )
                           flag_conv = 1; 
                       end
                        
                    end
                end
                
            end
        end
    end
end
        
end

