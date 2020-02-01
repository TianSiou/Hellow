function [ ew_new_1, ev_new, rsdl, cong_flag, ew_unit_c, ev_on_c, ew_inside_c_new, ev_inside_c_new, ...
           no_ew_in_c, next_ref_ev_outside_c, no_ref_ev, flag_conv, no_unit_ew ] = find_loss_real_ew( dim_n, ...
           mtx_NME, lambda_in, prev_ew, sigma, conv_radius, rad, rad_old, ew_new, ev_new, tolerance, ...
           ew_unit_c, no_ew_in_c, ew_inside_c_new, ev_on_c, ev_inside_c_new, next_ref_ev_outside_c, ...
           no_ref_ev, mtx_M, mtx_L, no_unit_ew )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ew_new_1    = [];
rsdl        = [];
cong_flag   = 0;
flag_conv   = 0;

tol_angle   = 1.0e-10;
rsdl_tol    = 1.0e-9;

idx_angle      = find(imag(lambda_in) >= -1.0e-10 & abs(lambda_in) <= 1);
tmp_ang2       = angle(lambda_in(idx_angle));
idxx           =  abs(tmp_ang2) <= tol_angle ;
tmp_ang2(idxx) = 0;
tmp_ang2       = mod(tmp_ang2,2*pi); 
 
t1           = angle(sigma);
if ( abs(t1) <= tol_angle )
    t1 = 0;
end

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
                            new_shift = abs(0.45*max_val_cur + 0.55*min_val_pre) * exp(1i*t1);
                        else
                            min_val_cur = min(real(real_neg_ew));
                            max_val_pre = max(real(prev_real_ew));
                            new_shift   = abs(0.55*min_val_cur + 0.45*max_val_pre) * exp(1i*t1);
                        end
                        
                        no_wanted_ew = 50;
                        
                        fprintf('\n Compute the eigenvalues in x-axis with shift = %11.4e+(%11.4e)*1i in find_loss_real_ew \n', real(new_shift), imag(new_shift));
                
                       [ deflation_enforce_1, lambda_in1, ew_new_1, ev_new, conv_radius_1, rsdl, no_deuplicate_ew_1, ...
                           cong_flag, ew_unit_c, ev_on_c, ew_inside_c_new, ev_inside_c_new, ...
                           no_ew_in_c, next_ref_ev_outside_c, no_ref_ev, idxx, no_unit_ew ] = compute_ew_near_sigma( dim_n, mtx_NME, ...
                           new_shift, 'no', rad, rad_old, ew_new, ev_new, no_wanted_ew, tolerance, ew_unit_c, ...
                           no_ew_in_c, ew_inside_c_new, ev_on_c, ev_inside_c_new, next_ref_ev_outside_c, ...
                           no_ref_ev, mtx_M, mtx_L, no_unit_ew );
                    
                       if ( isempty(rsdl) || max(rsdl) > rsdl_tol )
                           flag_conv = 1;
                       else
                           flag_conv = 0;
                       end
                        
                    end
                end
                
            end
        end
        
end

