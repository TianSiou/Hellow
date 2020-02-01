function [ dist_shift2ew, esti_conv_radius, new_shift, idx, dist_shift2_leftew ] = choose_shift_GTSHIRA( ...
    ew_S_Sinv, shift, idx, idx_ew, sigma, ew_skew_Cayley, ew_skew_Hamil, lambda_in, left_conv_point, ref_ew, ...
    prev_rad, pre_shift )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
            
            idx_target_lambda_in =  real(lambda_in) <= real(pre_shift) ;
            lambda_in            = lambda_in(idx_target_lambda_in); 
            
            tol_angle   = 1.0e-10;
            mu_0        = shift + 1/shift;
            angle_shift = mod(angle(shift), 2*pi);
            %
            % Compute the shortest distance between the convergent
            % eigenvalues in the previous radius and the shift value under the S+S^{-1}
            % transformation
            %
            delta_theta       = 20 / 180 * pi;
            angle_ref_ew      = angle(ref_ew);
            idx_tmp           = find(abs(angle_ref_ew) <= tol_angle);
            if ( ~isempty(idx_tmp) )
                angle_ref_ew(idx_tmp) = 0;
            end
            angle_ref_ew = mod(angle_ref_ew, 2*pi);
            
            %
            % Find the convergent eigenvalues which are near to the shift value
            %
            idx_2         = find(((angle_ref_ew >= angle_shift-delta_theta & angle_ref_ew <= angle_shift+delta_theta) ...
                | (real(ref_ew) <= real(abs(shift)*exp(1i*(angle_shift-delta_theta))) & ...
                real(ref_ew) >= real(abs(shift)*exp(1i*(angle_shift+delta_theta))))) & abs(ref_ew) <= prev_rad);
            if ( isempty(idx_2) )
                dist_shift2ew = 0;
            elseif (isempty(find(abs(angle_ref_ew - pi) > 1.5/180*pi, 1)))
                dist_shift2ew = 0;
            else
                %ref_ew(abs(angle_ref_ew - pi) > 1.5/180*pi)
                dist_shift2ew = min(abs(ew_S_Sinv(idx_2) - mu_0)); 
            end
            
            %
            % Compute the shortest distance between the convergent
            % eigenvalues in the previous iteration and the shift value under the S+S^{-1}
            % transformation
            %
            
            if ( isempty(lambda_in) )
                dist_shift2_leftew = 0;
            else                
                angle_lambda_in = angle(lambda_in);
                idx_tmp         = find(abs(angle_lambda_in) <= tol_angle);
                if ( ~isempty(idx_tmp) )
                    angle_lambda_in(idx_tmp) = 0;
                end
                angle_lambda_in = mod(angle_lambda_in,2*pi);
            
                idx_angle_lambda = find(abs(mod(angle(lambda_in),2*pi) - pi) > 1/180*pi, 1);
                %lambda_in(idx_angle_lambda)
                if ( isempty(idx_angle_lambda) || angle_shift - max(angle_lambda_in) > 50/180*pi)
                    dist_shift2_leftew = 0;
                else
                    tt2                = lambda_in + 1./lambda_in;
                    dist_shift2_leftew = min(abs(tt2 - mu_0));
                end
            end
            
            if ( ~isempty(left_conv_point) )
                dist_shift2_leftew = min(dist_shift2_leftew, abs(mu_0 - left_conv_point - 1 / left_conv_point));
            end
             
            %
            % Estimate the convergent radius of the GTSHIRA based on the
            % value produced by skew-J-Lanczos
            %
            % (i) remove the values which are too far away the shift value
            tt            = abs(ew_skew_Hamil(idx_ew(idx),1)-sigma);
            [~,idx_sort]  = sort(tt);
            leng_tt       = length(tt);
            tt2           = tt(idx_sort);
            tt3           = abs(tt2(2:leng_tt) - tt2(1:leng_tt-1));
            [dd,idx_max]  = max(tt3); 
            if (leng_tt >= 3 )
                if ( idx_max >= 2 )
                    if ( dd > 3 * max(tt3(1:idx_max-1))&& leng_tt - idx_max <= 2 )
                        fprintf('%2.0f estimated values are removed\n', leng_tt-idx_max);
                        idx = idx(idx_sort(1:idx_max));
                    end
                end
            end
            
            % (ii) Estimate the convergent radius esti_conv_radius 
            tt               = ew_skew_Cayley(idx,1) + 1./ew_skew_Cayley(idx,1);
            esti_conv_radius = max(abs(tt - mu_0));
            ratio            = 1.2;  
           
            % =========================
            % add new condition at 2016/12/9
            %
            if ( esti_conv_radius >= dist_shift2_leftew && dist_shift2ew / esti_conv_radius > 1.4 ) 
                dist_shift2ew = 0;
            end
            % 
            % ==========================
            
            idx_tt    = find(abs(ew_S_Sinv - mu_0) < ratio*esti_conv_radius);
            first_iew = length(idx_tt);
            %fprintf('Number of ew inside convergent radius is %2.0f \n',first_iew)
            fprintf('esti_conv_radius = %11.4e, dist_shift2ew = %11.4e, dist_shift2_leftew = %11.4e \n', ...
                esti_conv_radius, dist_shift2ew, dist_shift2_leftew);
            
%             if ( first_iew > 50 )
            if ( esti_conv_radius >= 0.98 * dist_shift2ew && esti_conv_radius >= 0.98 * dist_shift2_leftew)
                test_vec                = ew_S_Sinv(idx_2); %ew_S_Sinv(idx_tt,1);
                %no_test   = zeros(10,1);
                test_dist_shift2ew      = zeros(20,1);
                test_esti_conv_radius   = zeros(20,1);
                test_dist_shift2_leftew = zeros(20,1);
                shift_vec               = shift * linspace(0.95, 1, 20)'; 
                for i2 = 1:0 %length(shift_vec) 
                    mu_0                          = shift_vec(i2,1) + 1/shift_vec(i2,1);
                    if ( ~isempty(idx_2) )
                        test_dist_shift2ew(i2,1)  = min(abs(test_vec - mu_0));
                    end
                    test_esti_conv_radius(i2,1)   = max(abs(tt - mu_0));
                    if ( isempty(lambda_in) )
                        test_dist_shift2_leftew(i2,1) = abs(mu_0 - left_conv_point - 1 / left_conv_point);
                    else
                        test_dist_shift2_leftew(i2,1) = min(min(abs(tt2 - mu_0)), abs(mu_0 - left_conv_point - 1 / left_conv_point));
                    end
                    
%                     idd = find(abs(test_vec - mu_0) < ratio*esti_conv_radius);
%                     if ( ~isempty(idd) )
%                         no_test(i2,1)   = length(idd);
%                     end
                end
                
                idxx2       = find(test_esti_conv_radius > test_dist_shift2ew & test_esti_conv_radius > test_dist_shift2_leftew);
                if ( ~isempty(idxx2) )
                    new_shift          = shift_vec(idxx2(1));
                    dist_shift2ew      = test_dist_shift2ew(idxx2(1));
                    esti_conv_radius   = test_esti_conv_radius(idxx2(1));
                    dist_shift2_leftew = test_dist_shift2_leftew(idxx2(1));
                else
                    new_shift = shift;
                end
                
            else
                new_shift = shift;
            end
            
end

