function [ ew_unit_c, ev_on_c, ew_inside_c_new, ev_inside_c_new, no_deuplicate_ew, ...
    rad, ref_ev_outside_c, ref_ew_inside_c, ref_ev_inside_c, max_wanted_ew, flag_conv, ew_inner, ...
    radius_case, weighting_choosing, jj ] = compute_1st_circle_ew_0101( dim_n, mtx_NME, mtx_M, mtx_L, ...
    rad, rad_old, max_wanted_ew_org, radius_case, vidObj, weighting_choosing, ew_unit_c, ev_on_c, ...
    ew_inside_c, ref_ev_outside_c, ref_ew_inside_c, ref_ev_inside_c, plot_flag, deflation_flag )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
            
global ew 

ew_inside_c_new                 = zeros(dim_n,1);
ev_inside_c_new                 = zeros(2*dim_n,dim_n);
no_ew_in_c                      = length(ew_inside_c);
no_ew_in_c_org                  = no_ew_in_c;
% ew_inside_c_new(1:no_ew_in_c,1) = ew_inside_c;
% ev_inside_c_new(:,1:no_ew_in_c) = ev_inside_c;
% ref_ew_inside_c                 = ew_inside_c;
% ref_ev_inside_c                 = ev_inside_c; 
next_ref_ev_outside_c           = zeros(2*dim_n,floor(dim_n/5)); 
ew_S_Sinv                       = ew_inside_c + 1./ew_inside_c;
% no_ew                           = no_ew_in_c; 
ref_ew                          = ew_inside_c;
no_unit_ew                      = size(ew_unit_c,1);

idx_target_ref_ew   = find( abs(ref_ew) <= rad_old );
ref_ew              = ref_ew(idx_target_ref_ew);
ew_S_Sinv           = ew_S_Sinv(idx_target_ref_ew);

idx_neg_real_ew     = find(abs(imag(ew_inside_c)) <= 1.0e-8 & mod(angle(ew_inside_c),2*pi) >= pi - 3.0e-2 ...
    & mod(angle(ew_inside_c),2*pi) <= pi + 3.0e-2 & abs(ew_inside_c) >= rad);

% idx_neg_real_ew     = find(angle(ew_inside_c) >= pi - 3.0e-2 ...
%                 & angle(ew_inside_c) < pi + 1.0e-3 & abs(ew_inside_c) >= rad ); %rad_old - (rad_old - rad)/3);
if ( isempty(idx_neg_real_ew) )
    min_rad_neg_real_ew = -rad_old;
else
    min_rad_neg_real_ew = max(real(ew_inside_c(idx_neg_real_ew)));
end
   
max_jj              = 50;
leng_sigma          = 4;
max_iter_Bi_Lanczos = 200; 
NO_ew_T             = zeros(max_jj,1);
ew_skew_Hamil       = zeros(4*max_iter_Bi_Lanczos,leng_sigma);
ew_skew_Cayley      = zeros(4*max_iter_Bi_Lanczos,leng_sigma);
jj                  = 1; 

shift        = zeros(max_jj,1);
%rad          = 0.9 * rad;
% tmp_real     = rad * (1 - (1.0e-3));
% tmp_imag     = sqrt(rad^2 - tmp_real^2);
% shift(jj,1)  = tmp_real + 1i*tmp_imag;
idx_pos_ew_prev      = find( mod(angle(ref_ew_inside_c), 2*pi) < 3/180*pi );% & abs(ref_ew_inside_c) >= rad & abs(ref_ew_inside_c) <= rad_old );
leng_idx_pos_ew_prev = length(idx_pos_ew_prev);
 
% if ( leng_idx_pos_ew_prev == 0 )
%     shift(jj,1)  = rad * exp(1i * 3.5 / 180 * pi);
% else
if ( leng_idx_pos_ew_prev <= 20 && rad >= 1.0e-3)
    shift(jj,1)  = rad * exp(1i * 3 / 180 * pi);
elseif ( leng_idx_pos_ew_prev <= 30 && rad >= 1.0e-3 )
    shift(jj,1)  = rad * exp(1i * 1.5 / 180 * pi); %exp(1i * 3 / 180 * pi);
else
    shift(jj,1)  = rad * exp(1i * 1.0 / 180 * pi);
end
 
ew_inner      = zeros(max_jj, 1);
prediction_pt = zeros(max_jj,1);
shift_1       = zeros(max_jj,1);
shift_1(jj,1) = (shift(jj,1)+1)./(shift(jj,1)-1);
 
maxit         = 70;

shrink            = 'no'; 
arrive_right_tag  = 'no';
arrive_left_tag   = 'no';
tolerance         = 1.0e-12;

% idx_ew_in_c       = zeros(max_jj, 1);
% idx_ew_in_c(1,1)  = no_ew_in_c;
no_deuplicate_ew  = zeros(max_jj,2);
deflation_enforce = 'no';
deflate_tol       = 1.0e-9; %1.0e-10;
theta             = 0:0.001:2*pi;
x                 = cos(theta)+1i*sin(theta);
no_ref_ev         = 0; 
ratio             = 1.0; % less than or equal to 1
conv_radius       = zeros(max_jj,1);
% no_circle         = 1;
conv_region_angle = 0;
eigenwanted       = 10;
%maxit_GTSHIRA     = 100;

% figure(2) 
% plot(rad*x,'r-'); 
% drawnow
% hold on
  
if ( rad >= 1.0e-5 )
    max_no_re_estimate = 1;
    max_no_skew_J_Lan  = 2;
else
    max_no_re_estimate    = -1;
    max_no_skew_J_Lan     = 1;
    sigma                 = shift(jj,1);
end
%m2            = 8; %3
max_wanted_ew = max_wanted_ew_org;
prev_ew       = [];
tol_angle     = 1.0e-10;
no_ew_isempty = 0;
ve_ew_isempty = zeros(20,1);

if ( strcmp(weighting_choosing,'left') )
    checking_choice = 'yes';
else
    checking_choice = 'no';
end

flag_stop              = 'no';
ratio_rad              = 0.98;
prev_deflation_enforce = 'no';
rat_max_angle          = 1.55;

while ( jj <= max_jj && strcmp(shrink,'no') && strcmp(flag_stop,'no') ) %20
    skill_flag     = 0;
    no_re_estimate = 0;
    overlapping    = 'no';
    Move_shift     = 'no';
            
    while ( strcmp(overlapping,'no') && no_re_estimate <= max_no_re_estimate && strcmp(flag_stop,'no') )
        re_estimate_check = 'yes';
        no_skew_J_Lan     = 1;
    
        while ( no_skew_J_Lan <= max_no_skew_J_Lan && strcmp(re_estimate_check,'yes') && strcmp(flag_stop,'no') )
            sigma               = shift_1(jj, 1); %shift_1(jj, 1);
            %fprintf('\n skew_J_Lan with shift = %11.4e+(%11.4e)*1i and angle = %11.4e \n', real(shift(jj,1)), imag(shift(jj,1)), angle(shift(jj,1)));
            
            %[ LU_W_beta, ew_T ] = estimate_by_skew_J_Lanczos( sigma, dim_n, mtx_NME, maxit );
        
            shift(jj, 1)       = abs(shift(jj, 1)) * exp(1i * angle(shift(jj, 1))) * ratio_rad;
            sigma_old          = shift(jj, 1);
            fprintf('Run estimation I with sigma = %24.16e+(%24.16e)*1i \n', real(shift(jj, 1)), imag(shift(jj, 1)));
            W                  = mtx_NME.mtx_A - shift(jj, 1)*mtx_NME.mtx_Q + (shift(jj, 1)^2)*(mtx_NME.mtx_A.');
            LU_W.Perm_amd_vec  = amd(W);
            W_reorder          = W(LU_W.Perm_amd_vec, LU_W.Perm_amd_vec);
            LU_W.Perm_amd      = sparse(LU_W.Perm_amd_vec,1:dim_n,ones(dim_n,1));
        
            [LU_W.Low_L,LU_W.upper_U,LU_W.Perm_LU] = lu(W_reorder);
        
            LU_W.Perm_LU_Tran  = LU_W.Perm_LU.';
            LU_W.Low_L_Tran    = LU_W.Low_L.';
            LU_W.upper_U_Tran  = LU_W.upper_U.'; 
            LU_W.sigma         = shift(jj, 1);
        
            if ( abs(shift(jj, 1)) <= 3.0e-3 )
                LU_W.mtx_Qp  = W;
                LU_W.mtx_QpT = W.';
            end
            
            nmax               = min(150, maxit + round(0.9 / abs(shift(jj, 1))));

            [ ew_T, residual, R, H, Y, Z, nrm_Err_N ] = estimate_by_GTSHIRA( mtx_NME.mtx_A, mtx_NME.mtx_Q, ...
                shift(jj, 1), mtx_NME, @(x,transp_flag)solve_LS_PQEP(x,transp_flag,LU_W), dim_n, nmax );
            
            if ( jj == 1 )
                radius_estimate  = max(abs(ew_T + 1./ew_T - shift(jj, 1) - 1/shift(jj, 1)));
                dist_shift_neg_1 = abs(shift(jj, 1) + 1/shift(jj, 1) + 2); 
                if ( radius_estimate > dist_shift_neg_1 )
                    flag_stop = 'no'; %'yes';
                end
            end
            
            restart_info = 'yes';
            
            %fprintf('Number of estimated eigenvalues = %3.0f with %3.0f iterations\n',length(ew_T), nmax)
            
            if ( strcmp(plot_flag,'yes') )
                figure(4)
                plot(x,'r-')
                hold on
                plot(ew,'bo','LineWidth',1);
            
                plot(ew_T,'r+','LineWidth',2);
                plot(real(shift(jj,1)),imag(shift(jj,1)),'k^','LineWidth',2); 
                min_real = min(min(real(ew_T)), real(shift(jj,1)));
                max_real = max(max(real(ew_T)), real(shift(jj,1)));
                min_imag = min(min(imag(ew_T)), imag(shift(jj,1)));
                max_imag = max(max(imag(ew_T)), imag(shift(jj,1)));
                min_real = min_real * (1-sign(min_real)*0.1);
                max_real = max_real * (1+sign(max_real)*0.1);
                min_imag = min_imag * (1-sign(min_imag)*0.1);
                max_imag = max_imag * (1+sign(max_imag)*0.1);
                axis([min_real max_real min_imag max_imag])
                hold off
                drawnow
            end
            
            angle_shift  = angle(shift(jj,1));
            angle_1st_ew = angle(ew_T(1,1));
            if ( abs(angle_shift) < 1.0e-7 )
                angle_shift = 0;
            end
            if ( abs(angle_1st_ew) < 1.0e-7 )
                angle_1st_ew = 0;
            end
            
            %if ( (abs(ew_T(1,1)) >= rad_old || (angle_1st_ew - angle_shift) > 50/180*pi ) && residual(1) < 1.0e-2)
            if ( (angle_1st_ew - angle_shift) > 50/180*pi  && residual(1) < 1.0e-2)
                cn_radius   = abs(ew_T(1,1) + 1 / ew_T(1,1) - shift(jj,1) - 1 / shift(jj,1));
                idx_pos_img = find( imag(ew_T) > 0 );
                ang_thm_pt  = esstimating_fm_in_conv_circle(cn_radius, shift(jj,1), ew_T(idx_pos_img(1)), 'left');
                ang_thm_pt  = min( 179/180*pi, ang_thm_pt);
                shift(jj,1) = rad * exp(1i * ang_thm_pt); 
                sigma       = shift(jj,1);
                Move_shift  = 'yes';
                
                fprintf('Moving shift value to %11.4e + (%11.4e) * 1i \n', real(shift(jj,1)), imag(shift(jj,1)));
                if ( abs(angle_shift - 179/180*pi) < 1.0e-7 )
                    %re_estimate_check = 'no';
                    shrink            = 'yes';
%                 else
% %                     no_skew_J_Lan = max_no_skew_J_Lan;
%                 end
                else
                    nrm_Err_N     = 1;
                    ij            = 1;
                    tol_nrm_Err_N = 1.0e-3;
                    while ( ij <= 2 && nrm_Err_N > tol_nrm_Err_N )
                        re_estimate_check = 'no';
                        
                        shift(jj, 1)      = abs(shift(jj, 1)) * exp(1i * angle(shift(jj, 1))) * ratio_rad;
                        sigma_old         = shift(jj, 1);
                
                        fprintf('Run re_estimation II with sigma = %24.16e+(%24.16e)*1i \n', real(shift(jj, 1)), imag(shift(jj, 1)));
                        W                  = mtx_NME.mtx_A - shift(jj, 1)*mtx_NME.mtx_Q + (shift(jj, 1)^2)*(mtx_NME.mtx_A.');
                        LU_W.Perm_amd_vec  = amd(W);
                        W_reorder          = W(LU_W.Perm_amd_vec, LU_W.Perm_amd_vec);
                        LU_W.Perm_amd      = sparse(LU_W.Perm_amd_vec,1:dim_n,ones(dim_n,1));
        
                        [LU_W.Low_L,LU_W.upper_U,LU_W.Perm_LU] = lu(W_reorder);
        
                        LU_W.Perm_LU_Tran  = LU_W.Perm_LU.';
                        LU_W.Low_L_Tran    = LU_W.Low_L.';
                        LU_W.upper_U_Tran  = LU_W.upper_U.'; 
                        LU_W.sigma         = shift(jj, 1);
        
                        if ( abs(shift(jj, 1)) <= 3.0e-3 )
                            LU_W.mtx_Qp  = W;
                            LU_W.mtx_QpT = W.';
                        end
            
                        nmax               = min(150, maxit + round(0.9 / abs(shift(jj, 1))));

                        [ ew_T, residual, R, H, Y, Z, nrm_Err_N ] = estimate_by_GTSHIRA( mtx_NME.mtx_A, ...
                            mtx_NME.mtx_Q, shift(jj, 1), mtx_NME, ...
                            @(x,transp_flag)solve_LS_PQEP(x,transp_flag,LU_W), dim_n, nmax );
                    
                        ij = ij + 1;
                        if ( nrm_Err_N >= tol_nrm_Err_N )
                            shift(jj, 1) = rad * exp(1i*mod(angle(shift(jj, 1)), 2*pi)*1.005);
                        end
                    end
            
                    mm           = length(ew_T); 
                    eigenwanted  = mm; % - 10;
                    restart_info = 'yes';
                    sigma        = shift(jj,1);
                    
                end
                mm                = 0;
        
            else
                if ( ~isempty(ew_T) || jj == 1)
                    re_estimate_check = 'no';
                elseif ( jj > 1 )
                    idx2               = imag(lambda_in) >= 0;
                    right_angle        = max(angle(lambda_in(idx2)));
                    left_angle         = angle(shift(jj,1));
                    %[right_angle conv_region_angle left_angle left_angle - conv_region_angle 5/180*pi ]
                    if (left_angle >= pi - 3.0e-2 && abs(left_angle - right_angle) <= 10/180*pi ...
                            && strcmp(arrive_left_tag,'no') && strcmp(arrive_right_tag,'no') ) 
                        gamma_real        = -rad * 0.99; %-rad + 1e-3;
                        gamma_imag        = sqrt(rad^2 - gamma_real^2);
                        shift(jj,1)       = gamma_real + 1i * gamma_imag;
                        shift_1(jj, 1)    = (shift(jj, 1)+1)/(shift(jj, 1)-1); 
                        re_estimate_check = 'no';
                    elseif (abs(left_angle - conv_region_angle) >= 5/180*pi)  % original version right_angle
                        shift_angle       = left_angle - abs(left_angle - conv_region_angle) / 5;
                        shift(jj,1)       = rad * exp(1i*shift_angle);
                        shift_1(jj, 1)    = (shift(jj, 1)+1)/(shift(jj, 1)-1); 
                        fprintf('No estimated eigenvalues, give new angle %11.4e \n', shift_angle);
                        no_skew_J_Lan     = no_skew_J_Lan + 1;
                        %continue; conv_region_angle
                    else
                        re_estimate_check = 'no';
                    end
                    fprintf('The angle of shift value is %11.4e \n', left_angle);
                else
                    gamma_imag     = imag(shift(jj,1)) * 1.05;
                    gamma_real     = sqrt(rad^2 - gamma_imag^2);
                    shift(jj,1)    = gamma_real + 1i * gamma_imag;
                    shift_1(jj, 1) = (shift(jj, 1)+1)/(shift(jj, 1)-1);
                    no_skew_J_Lan  = no_skew_J_Lan + 1;
                
                    disp(shift(jj,1))
                    %re_estimate_check = 'no';
                end
                
                mm          = length(ew_T); 
                eigenwanted = mm; % - 10;
        
            end
        end
    
%         mm          = length(ew_T); 
%         eigenwanted = mm - 10;
        fprintf('number eigenvalues = %3.0f, iteration = %3.0f and sigma = %11.4e + (%11.4e)*1i \n', mm, nmax, real(shift(jj,1)), imag(shift(jj,1)));
        overlapping    = 'yes';
         
        if ( rad < 1.0e-3 && jj == 1 )
            sigma = shift(jj,1);
        elseif ( mm ~= 0 )
        
            ew_skew_Hamil(1:mm,jj)              = ew_T;
            idx_ew                              = find(abs(ew_T) < rad_old);
            NO_ew_T(jj, 1)                      = length(idx_ew);
            ew_skew_Cayley(1:NO_ew_T(jj, 1),jj) = ew_T(idx_ew, 1);
                    
            idx = find(imag(ew_skew_Cayley(1:NO_ew_T(jj, 1),jj)) >= 0);
        
            if (~isempty(idx))
                left_conv_point = rad * exp(1i*conv_region_angle);
                %
                % Determine the shift value for GTSHIRA
                %
                if ( jj == 1 )
                    [ dist_shift2ew, esti_conv_radius, new_shift, idx, dist_shift2_leftew ] = ...
                        choose_shift_GTSHIRA( ew_S_Sinv, shift(jj,1), idx, idx_ew, sigma, ew_T, ... %ew_skew_Cayley(:,jj), ...
                        ew_skew_Hamil(:,jj), [], [], ref_ew, rad_old, 0 ); %rad, rad, ref_ew, rad_old, 0 );
                else
                    [ dist_shift2ew, esti_conv_radius, new_shift, idx, dist_shift2_leftew ] = ...
                        choose_shift_GTSHIRA( ew_S_Sinv, shift(jj,1), idx, idx_ew, sigma, ew_T, ... %ew_skew_Cayley(:,jj), ...
                        ew_skew_Hamil(:,jj), new_left_ew, left_conv_point, ref_ew, rad_old, shift(jj-1,1) );
                end
            
                idxx_unit_ew        = abs(abs(ew_skew_Cayley(idx,jj))-1) < 1.0e-5;
                min_ang_skew_Cay_ew = min(angle(ew_skew_Cayley(idx(idxx_unit_ew),jj)));
                if ( jj == 1 )
                    max_ang_ew = 0;
                else
                    ang_tmp    = angle(new_left_ew);
                    idx_tmp2   = find( abs(ang_tmp) <= tol_angle );
                    if( ~isempty(idx_tmp2) )
                        ang_tmp(idx_tmp2) = 0;
                    end
                    ang_tmp    = mod(ang_tmp, 2*pi);
                    
                    max_ang_ew = max(conv_region_angle,max(ang_tmp));
                end
            
                %
                % If esti_conv_radius is less than dist_shift2ew, then maybe some of
                % eigenvalues will loss. It needs to re-estimate the shift
                % value 
                %
                re_estimate_flag = 'no';
                switch radius_case
                    case 'normal'
%                         fprintf('normal case \n');
%                         [ esti_conv_radius dist_shift2ew dist_shift2_leftew ratio ]
                        if (esti_conv_radius < ratio * dist_shift2ew || esti_conv_radius < ratio * dist_shift2_leftew)
                            re_estimate_flag = 'yes';
                        end
                        
                    case 'weighting'
%                         fprintf('weighting case \n');
%                         [ esti_conv_radius dist_shift2_leftew ratio ]
                        if (esti_conv_radius < ratio * dist_shift2_leftew)
                            re_estimate_flag = 'yes';
                        end
                        
                end
                
                if (strcmp(arrive_left_tag,'no') && strcmp(arrive_right_tag,'no') ...
                    && strcmp(re_estimate_flag,'yes') ) 
                
                    if ( ~(esti_conv_radius >= ratio * dist_shift2_leftew && dist_shift2ew > esti_conv_radius*10 ) ) 
                    % 
                    % Determine new shift value
                    %
                        diff_angle     = abs(min_ang_skew_Cay_ew - max_ang_ew);
                        ang_shift      = angle(shift(jj,1));
                        ang_shift2ew   = abs(ang_shift - max_ang_ew);
                    
                        if ( max_ang_ew >= ang_shift - 1.0e-3 && max_ang_ew >= pi - 3.0e-3 ) % spectial case
                            shift(jj,1) = -rad * 0.98 + 1i * 1.0e-3 * rad;
                        else
                            if ( diff_angle >= ang_shift2ew * 0.2 )
                                if ( jj == 1 )  
                                    new_theta = max(min(ang_shift - diff_angle / 2, min_ang_skew_Cay_ew), 1.0e-4);
                                else
                                    new_theta = max(min(pi,max_ang_ew+5/180*pi), min(ang_shift - diff_angle / 2, min_ang_skew_Cay_ew));
                                end
                                fprintf('new_theta = %11.4e, min_ang_skew_Cay_ew = %11.4e \n', new_theta, min_ang_skew_Cay_ew)
                            else
                                if ( jj == 1 )
                                    new_theta = min(ang_shift - diff_angle / 2, min_ang_skew_Cay_ew);
                                else 
                                    if ( max_ang_ew+5/180*pi < ang_shift)
                                        new_theta = max(min(pi,max_ang_ew+5/180*pi), ang_shift - diff_angle / 2);
                                    else
                                        new_theta = ang_shift - diff_angle / 2; 
                                    end 
                                end
                            end 
                            shift(jj,1)    = rad * exp(1i*new_theta); %gamma_real + 1i * gamma_imag;
                        end
                        shift_1(jj, 1) = (shift(jj, 1)+1)/(shift(jj, 1)-1); 
                        new_shift      = shift(jj,1);
                     
                    
                        if ( (ang_shift > max_ang_ew+1.0e-1 && no_re_estimate <= max_no_re_estimate) || (jj == 1 && rad >= 1.0e-3) )   
                            no_re_estimate        = no_re_estimate + 1;
                            fprintf('No overlapping, re-estimate with theta = %11.4e and rad = %11.4e!!\n', new_theta, rad);
                            overlapping           = 'no';
                            sigma                 = new_shift; 
                        end
                    end
                
                end
        
                if ( strcmp(overlapping,'yes') || no_re_estimate > max_no_re_estimate )
%                     figure(3)
%                     tt = [sigma, -sigma, conj(sigma), -conj(sigma)]; 
%                     plot(tt,'g+','LineWidth',2);
%                     figure(3)
%                     if ( ~isempty(idx_ew) && ~isempty(idx) )
%                         plot(ew_skew_Hamil(idx_ew(idx),jj),'rx','LineWidth',2);
%                     end
  
                    sigma                 = new_shift; %shift(jj, 1); 
                    mu_0                  = sigma + 1 / sigma;
      
                end
            else
                sigma                 = shift(jj,1); 
            end 
        
        elseif ( strcmp(arrive_left_tag,'no') && strcmp(arrive_right_tag,'no') && strcmp(Move_shift, 'no') )  
            %
            % Re-estimate new shift value 
            %
            if ( jj == 1 && rad <= 5.0e-2 )
                overlapping    = 'yes';
                sigma          = shift(jj,1);
            else
                if ( jj == 1 )
                    tmp_real       = rad - (1.0e-3);
                    tmp_imag       = sqrt(rad^2 - tmp_real^2);
                    shift(jj,1)    = tmp_real + 1i*tmp_imag; 
                else 

                    shift_angle       = angle(shift(jj,1));
                    shift(jj,1)       = rad * exp(1i * (shift_angle + abs(shift_angle - conv_region_angle) / 2));
                    conv_region_angle = shift_angle; 
                end
                shift_1(jj, 1) = (shift(jj, 1)+1)/(shift(jj, 1)-1);
                overlapping    = 'no';
                no_re_estimate = no_re_estimate + 1;
                sigma          = shift(jj,1);
                fprintf('Re-doing again to estimate approximate eigenvalues \n') %with length of target_ew_skew_Cayley %3.0f \n',length(target_ew_skew_Cayley));
                disp(shift(jj,1)) 
            end
        end
        %end
    end
   
    %
    % Based on no_ew_isempty, define the maximual iteration numbers of GTSHIRA
    %
    if ( no_ew_isempty <= 3 )
        maxit_GTSHIRA = 25; %100;
    else
        maxit_GTSHIRA = 30; %50;
    end
%
%  ============ Begin to compute eigenvalues near to the shift value ============
%
    if ( strcmp(flag_stop,'no') )

        if ( abs(sigma_old - shift(jj,1)) > 1.0e-12 )  
            sigma              = ratio_rad * abs(sigma) * exp(1i*angle(sigma));
            fprintf('Compute LU factroization with sigma = %24.16e+(%24.16e)*1i \n\n', real(sigma), imag(sigma));
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
            
            restart_info = 'no';
        end
           
        % ===============
        if ( (strcmp(arrive_left_tag,'yes') || strcmp(arrive_right_tag,'yes')) && rad < 1.0e-3 )
            nmax                = min(150, maxit + round(0.9 / abs(shift(jj, 1))));
            max_it              = 2;
            iter                = 1;
            re_estimate_real_ew = 'yes';

            while ( iter <= max_it && strcmp(re_estimate_real_ew, 'yes') )
                
                fprintf('Run real estimation with sigma = %24.16e+(%24.16e)*1i \n\n', real(sigma), imag(sigma));
                
                [ ew_T, residual, R, H, Y, Z, nrm_Err_N ] = estimate_by_GTSHIRA( mtx_NME.mtx_A, mtx_NME.mtx_Q, sigma, ...
                        mtx_NME, @(x,transp_flag)solve_LS_PQEP(x,transp_flag,LU_W), dim_n, nmax );
                    
                min_real_part = min(real(ew_T));
                max_real_part = max(real(ew_T));
                restart_info  = 'yes';
                
                if ( strcmp(plot_flag,'yes') )
                    figure(4)
                    plot(x,'r-')
                    hold on
                    plot(ew,'bo','LineWidth',1);
            
                    plot(ew_T,'r+','LineWidth',2);
                    plot(real(shift(jj,1)),imag(shift(jj,1)),'k^','LineWidth',2); 
                    min_real      = min(min_real_part, real(shift(jj,1)));
                    max_real      = max(max_real_part, real(shift(jj,1)));
                    min_imag      = min(min(imag(ew_T)), imag(shift(jj,1)));
                    max_imag      = max(max(imag(ew_T)), imag(shift(jj,1)));
                    min_real      = min_real * (1-sign(min_real)*0.1);
                    max_real      = max_real * (1+sign(max_real)*0.1);
                    min_imag      = min_imag * (1-sign(min_imag)*0.1);
                    max_imag      = max_imag * (1+sign(max_imag)*0.1);
                    axis([min_real max_real min_imag max_imag])
                    hold off
                    drawnow
                end
            
                %
                % Have any intersection for the estimating eigenvalues and the
                % convergent eigenvalues in the previous step
                % If yes, then do GTSHIRA
                % otherwise, re-estimate the shift value
                %
                estimate_conv_rad  = max(abs( ew_T + 1./ew_T - sigma - 1 / sigma));
                idx2               = imag(lambda_in) >= -1.0e-10;
                [min_val, idx_min] = min(abs(lambda_in(idx2) + 1./lambda_in(idx2) - sigma - 1 / sigma)); 
                leng_spectrum      = abs(max_real_part - min_real_part); 
                if ( min_val > estimate_conv_rad )
                    if ( strcmp(arrive_left_tag,'yes') )
                        leng_loss_interval = abs(real(lambda_in(idx2(idx_min))) - max_real_part);
                        if ( leng_loss_interval >= 0.85 * leng_spectrum )
                            tt = min(0.75, leng_spectrum / leng_loss_interval / 2);
                            tt = ( 1 - tt ) * max_real_part + tt * real(lambda_in(idx2(idx_min)));
                            %fprintf('case 1\n');
                        else
                            tt = 0.65;
                            tt = tt * real(lambda_in(idx2(idx_min))) + ( 1 - tt ) * real(sigma);
                            %fprintf('case 2\n');
                        end
                    else
                        leng_loss_interval = abs(real(lambda_in(idx2(idx_min))) - min_real_part);
                        if ( leng_loss_interval >= 0.85 * leng_spectrum )
                            tt = min(0.75, leng_spectrum / leng_loss_interval / 2);
                            tt = tt * min_real_part + (1 - tt) * real(lambda_in(idx2(idx_min)));
                            %fprintf('case 3\n');
                        else
                            tt = 0.65;
                            tt = tt * real(lambda_in(idx2(idx_min))) + ( 1 - tt ) * real(sigma);
                            %fprintf('case 4\n');
                        end
                    end
                    
                    if ( abs(abs(tt)-1) < 1.0e-9 )
                        sigma       = abs(tt) * exp(1i * angle(sigma)) * ratio_rad;
                    else
                        sigma       = abs(tt) * exp(1i * angle(sigma));
                    end
                    shift(jj,1) = sigma;
                    
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
                    
                else
                    re_estimate_real_ew = 'no';
                end
            
                iter = iter + 1;
            end
            
        end
        % ==============
         
        %  ======================
        % added at 2016/12/5
        %
        if ( strcmp(restart_info,'yes') )
            mu_0                 = sigma + 1/sigma;
            estimate_conv_radius = max(abs(ew_T + 1./ew_T -mu_0));
            dist_mu0_pos1          = abs(mu_0 - 2);
            dist_mu0_neq1          = abs(mu_0 + 2);
            if ( estimate_conv_radius >= dist_mu0_pos1 && estimate_conv_radius >= dist_mu0_neq1 )
                deflation_enforce = 'no';
            end
        end
        %  ======================
        %
        if ( strcmp(deflation_flag, 'no') )
            deflation_enforce = 'no';
        end
        
%         if ( strcmp(flag_stop,'no') )
            %if (jj > 1 && strcmp(deflation_enforce, 'yes') && real(sigma) > 0 )
            if (jj > 1 && strcmp(deflation_enforce, 'yes') && rad > 2.0e-3 ) % 6.0e-4
                
                % 
                % Find the eigenvalues to be in the convergent region
                % 

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
                    
                    if ( rad > 6.0e-3 )
                        no_defl_ew    = 20;
                    else
                        no_defl_ew    = 10;
                    end
                    
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
                
                fprintf('GTSHIRA with sigma = %11.4e + (%11.4e)*1i \n', real(sigma), imag(sigma));
                if ( flag == 0 )
                    restart_info  = 'no';
                    if ( real(sigma) >= -0.01 )
                        [ ew_new, ev_new, outer_it ] = Drive_GTSHIRA( dim_n, ...
                             @(x,transp_flag)Defl_mtx_vec_A( x, transp_flag, mtx_NME, mtx_deflate ), ...
                             @(x)Defl_mtx_vec_Q( x, mtx_NME, mtx_deflate ), mtx_NME, ...
                             @(x,transp_flag)solve_LS_PQEP(x,transp_flag,LU_W, mtx_deflate), ...
                             sigma, eigenwanted, tolerance, maxit_GTSHIRA, restart_info, R, H, Y, Z );
                    else
                        
                        [ ew_new, ev_new, outer_it ] = Drive_GTSHIRA( dim_n, ...
                             @(x,transp_flag)Defl_mtx_vec_A( x, transp_flag, mtx_NME, mtx_deflate ), ...
                             @(x)Defl_mtx_vec_Q_pone( x, mtx_NME, mtx_deflate ), mtx_NME, ...
                             @(x,transp_flag)solve_LS_PQEP(x,transp_flag,LU_W, mtx_deflate), ...
                             sigma, eigenwanted, tolerance, maxit_GTSHIRA, restart_info, R, H, Y, Z );
                    
                    end
                    
                    prev_deflation_enforce = 'yes';
                    deflation_enforce      = 'no';
                else
                    fprintf('GTSHIRA with sigma = %11.4e + (%11.4e)*1i \n', real(sigma), imag(sigma));
                    [ ew_new, ev_new, outer_it ] = Drive_GTSHIRA( dim_n, mtx_NME.mtx_A, mtx_NME.mtx_Q, mtx_NME, ...
                        @(x,transp_flag)solve_LS_PQEP(x,transp_flag,LU_W), sigma, eigenwanted, tolerance, ...
                        maxit_GTSHIRA, restart_info, R, H, Y, Z );
                    
                    deflation_enforce      = 'yes';
                    prev_deflation_enforce = 'no';
                end

            else 
                fprintf('GTSHIRA with sigma = %11.4e + (%11.4e)*1i \n', real(sigma), imag(sigma));
                [ ew_new, ev_new, outer_it ] = Drive_GTSHIRA( dim_n, mtx_NME.mtx_A, mtx_NME.mtx_Q, mtx_NME, ...
                    @(x,transp_flag)solve_LS_PQEP(x,transp_flag,LU_W), sigma, eigenwanted, tolerance, ...
                    maxit_GTSHIRA, restart_info, R, H, Y, Z );
                
                if ( strcmp(deflation_enforce, 'no') )
                    deflation_enforce      = 'yes';
                    prev_deflation_enforce = 'no';
                end

            end
             
            
            %fprintf('Before, no_ew_in_c = %4.0f \n', no_ew_in_c); 
            if ( ~isempty(ew_new) ) % GTSHIRA converges or not
                
                %
                flag_cong_GTSHIRA = 0;
                mu_0              = sigma + 1 / sigma;
                conv_radius(jj,1) = max(abs(ew_new + 1./ew_new - mu_0));
                dist_shift2neg1   = abs(mu_0 + 2);
                fprintf('conv_radius = %11.4e, dist_shift2neg1 = %11.4e \n',conv_radius(jj,1), dist_shift2neg1); 
                if ( conv_radius(jj,1) > dist_shift2neg1)
                    shrink     = 'yes';
                end
                    
                %
                % Find the deuplicated eigenvalues
                %
                leng_ew       = length(ew_new);
                lambda_in     = ew_new(1:2:leng_ew); 
                
                if ( ~isempty(find(abs(lambda_in) > 1+1.0e-8, 1)) )
                    fprintf('Error in choosing lambda_in \n');
                    flag_conv = 1;
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
                    if ( ~isempty(ew_unit_c) )
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
                    else
                        idx_new_ui_ew(idx_unit,1) = 1;
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
                        if (min(abs(lambda_in(idx_inside(i1))-ew_inside_c_new(interval,1))) < tol )
                            idx_deupli(idx_inside(i1),1) = 1;
                        end
                    end
                
                end
            
                idxx = find(idx_new_ui_ew == 1);
                if ( ~isempty(idxx) )
                    mm3                               = length(idxx); 
                    %tt                                = lambda_in(idxx)
                    %size(tt)
                    ew_unit_c(no_unit_ew+1:no_unit_ew+mm3,1) = lambda_in(idxx);
                    ew_unit_c(no_unit_ew+1:no_unit_ew+mm3,2) = ew_new(2*idxx);
                    %ew_unit_c(:,1)                    = [ew_unit_c(:,1); lambda_in(idxx)];
                    %ew_unit_c(:,2)                    = [ew_unit_c(:,2); ew_new(idxx+1)];
                    ev_on_c(:,no_unit_ew+1:no_unit_ew+mm3,1) = ev_in_c(:,idxx);
                    ev_on_c(:,no_unit_ew+1:no_unit_ew+mm3,2) = ev_out_c(:,idxx);
                    no_unit_ew                               = no_unit_ew + mm3;
                end
            
                idxx                   = find(idx_deupli == 0 & imag(lambda_in) > -1.0e-10); 
                mm2                    = length(idxx);
                percentage             = (leng_ew/2-mm2)/(leng_ew/2);
                no_deuplicate_ew(jj,:) = [ leng_ew/2-mm2  (leng_ew/2-mm2)/(leng_ew/2)];
                fprintf('%2.0f-th, percentage of deuplicated eigenvalues = %11.4e \n',jj, percentage)
                
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
                
                %
                % Counting how many times (continuous) that all computed eigenvalues are
                % duplicated.
                %
                prev_no_ew_isempty = no_ew_isempty;
                if ( mm2 == 0 )
                    no_ew_isempty       = no_ew_isempty + 1;
                else
                    ve_ew_isempty(jj,1) = 1;
                    no_ew_isempty       = 0;
                end
                
                no_ew_new = length(ew_new);
                rsdl      = ones(no_ew_new/2,1);
                for ki = 1:2:no_ew_new
                    if abs(imag(ew_new(ki)))<1e-10
                        ew_new(ki, 1) = real(ew_new(ki, 1));
                        ev_new(:, ki) = real(ev_new(:, ki));
                    end
                    ev_new(:,ki)     = ev_new(:,ki) / norm(ev_new(:,ki));
                    rsdl_vec         = mtx_M * ev_new(:,ki) - ew_new(ki) * (mtx_L * ev_new(:,ki)); 
                    rsdl((ki+1)/2,1) = norm(rsdl_vec);
                    fprintf('rsdl(%2.0f,%24.16e+(%24.16e)1i) = %10.4e \n', (ki+1)/2, real(ew_new(ki)), ...
                    imag(ew_new(ki)), rsdl((ki+1)/2,1));
%                     if ( rsdl > 1.0e-7 )
%                         skill_flag = 1;
% %                         flag_conv = 1;
% %                         return
%                     end
                end
                
                idx_rsdl          = find( rsdl <= 1.0e-9 );
                leng_idx_rsdl     = length(idx_rsdl);
                
                if ( rad >= 0.01 && leng_idx_rsdl < no_ew_new/2 )
                    flag_redo_GTSHIRA = 'yes';
                    if ( strcmp(prev_deflation_enforce,'no') )
                        sigma             = rad * exp(1i * mod(angle(sigma), 2*pi) * 0.993);
                    else
                        deflation_enforce = 'no';
                        no_unit_ew        = no_unit_ew - mm3;
                    end
                else
                    flag_redo_GTSHIRA = 'no';
                end
                    
                if ( leng_idx_rsdl == 0 )
                    skill_flag = 1;
                    if ( angle(sigma) >= 176/180*pi )
                        shrink     = 'yes';
                        fprintf('shrink in compute_jth_circle_ew \n');
                    end
                    
                elseif ( leng_idx_rsdl ~= mm2 )
                    
                    mm2               = leng_idx_rsdl;
                    idxx              = idxx(idx_rsdl); 
                    %
                    idx2              = zeros(mm2*2,1); 
                    idx2(1:2:2*mm2,1) = 2*idx_rsdl-1;
                    idx2(2:2:2*mm2,1) = 2*idx_rsdl; 
                    ew_new            = ew_new(idx2);
                    ev_new            = ev_new(:,idx2); 
                    
                end
                
                if ( skill_flag == 0 && strcmp(flag_redo_GTSHIRA,'no') )
                    ew_inside_c_new(no_ew_in_c+1:no_ew_in_c+mm2,1)     = lambda_in(idxx);
                    ev_inside_c_new(:,no_ew_in_c+1:no_ew_in_c+mm2)     = ev_in_c(:,idxx);
                    no_ew_in_c                                         = no_ew_in_c + mm2;  
                    next_ref_ev_outside_c(:,no_ref_ev+1:no_ref_ev+mm2) = ev_out_c(:,idxx);
                    no_ref_ev                                          = no_ref_ev + mm2;
                elseif (strcmp(flag_redo_GTSHIRA,'yes'))
                    
                    deflation_enforce = 'no';
                    
                    %fprintf('no_ew_in_c = %4.0f \n', no_ew_in_c);
                    fprintf('\n Re_doing GTSHIRA with sigma = %11.4e +(%11.4e)*1i \n', real(sigma), imag(sigma));
                    [ deflation_enforce, lambda_in, ew_new, ev_new, conv_radius(jj,1), rsdl, ...
                            no_deuplicate_ew(jj,:), cong_flag, ew_unit_c, ev_on_c, ew_inside_c_new, ...
                            ev_inside_c_new, no_ew_in_c, next_ref_ev_outside_c, no_ref_ev, idxx, no_unit_ew ] = ...
                            compute_ew_near_sigma( dim_n, mtx_NME, sigma, deflation_enforce, ...
                            rad, rad_old, ew_new, ev_new, eigenwanted, tolerance, ew_unit_c, ...
                            no_ew_in_c, ew_inside_c_new, ev_on_c, ev_inside_c_new, ...
                            next_ref_ev_outside_c, no_ref_ev, mtx_M, mtx_L, no_unit_ew );
                        
                    deflation_enforce = 'no';
                end
                
            else
                prev_no_ew_isempty = no_ew_isempty;
                flag_cong_GTSHIRA  = 1;
                no_ew_isempty      = no_ew_isempty + 1;
                
                if ( radius_estimate > dist_shift_neg_1 )
                    flag_stop = 'yes';
                end
            end
%
%  ============ End of computing eigenvalues near to the shift value ============
% 
            if ( ~isempty(ew_new) )
                ew_no_dupli = lambda_in(idxx);
                %if ( isempty(ew_no_dupli) )
                    checking_angle = mod(angle(sigma),2*pi);
                %else
                %    checking_angle = min(mod(angle(ew_no_dupli),2*pi));
                %end
             
                if ( flag_cong_GTSHIRA == 0 && checking_angle < 179/180*pi && ~isempty(ew_no_dupli) && skill_flag == 0 )
                %if ( flag_cong_GTSHIRA == 0 && mod(angle(sigma),2*pi) < 179/180*pi )
                    if ( rad < 1.0e-3 || strcmp(deflation_flag, 'no') )
                        deflation_enforce_new = 'no';
                    else
                        deflation_enforce_new = deflation_enforce;
                    end
            
                    [ ew_new_1, ev_new_1, rsdl, cong_flag, ew_unit_c, ev_on_c, ew_inside_c_new, ev_inside_c_new, ...
                        no_ew_in_c, next_ref_ev_outside_c, no_ref_ev, flag_conv, arrive_left_tag, no_unit_ew ] = ...
                        checking_computing_loss_ew( dim_n, mtx_NME, ew_no_dupli, prev_ew, ew_new, ... %lambda_in, prev_ew, ew_new, ...
                        ev_new, ref_ew_inside_c, sigma, eigenwanted, tolerance, conv_radius(jj,1), rad, ...
                        rad_old, ew_unit_c, no_ew_in_c, ew_inside_c_new, ev_on_c, ev_inside_c_new, ...
                        next_ref_ev_outside_c, no_ref_ev, mtx_M, mtx_L, ew_inside_c_new(no_ew_in_c-mm2+1:no_ew_in_c,1), ...
                        deflation_enforce_new, min_rad_neg_real_ew, arrive_left_tag, no_unit_ew);
                
                else
                    ew_new_1 = [];
                end
            
                prev_ew = lambda_in(idxx); %lambda_in;
            end
            %
            % Choose next shift value
            %
            if ( flag_cong_GTSHIRA == 0 )
                    tol_imag        = -1.0e-10; %0;
                    
                    idx_neg_real_ew = find(abs(imag(ew_inside_c_new(no_ew_in_c_org+1:no_ew_in_c,1))) <= 1.0e-5 & ...
                        real(ew_inside_c_new(no_ew_in_c_org+1:no_ew_in_c,1)) <= -1.0e-7);
                    
                    idx_ew          = find(imag(ew_new) >= tol_imag & abs(ew_new) <= 1); 
                    if ( isempty(idx_ew) )
                        fprintf('All eigenvalues are deuplicated \n');
                                                                        
                        idxx2       = find(imag(lambda_in) >= -1.0e-11);
                        increment   = 'no';
                        
                        % Let
                        %      C = {x; | x - sigma - 1/sigma | <= conv_radius}
                        % be the convergent region of GTSHIRA with shift sigma.
                        %
                        % Find the angle 'conv_region_angle' of x = rad * exp(1i*conv_region_angle) satisfying
                        %               | x - sigma - 1/sigma | = conv_radius
                        % 
%                         tmp         = [ (-(conv_radius(jj,1)-mu_0)+sqrt((conv_radius(jj,1)-mu_0)^2-4))/(2*rad); 
%                              (-(conv_radius(jj,1)-mu_0)-sqrt((conv_radius(jj,1)-mu_0)^2-4))/(2*rad);
%                              ((conv_radius(jj,1)+mu_0)+sqrt((conv_radius(jj,1)+mu_0)^2-4))/(2*rad);
%                              ((conv_radius(jj,1)+mu_0)-sqrt((conv_radius(jj,1)+mu_0)^2-4))/(2*rad) ];
%                         idx_tmp     = imag(tmp) >= -1.0e-10; 
                        ang_thm_pt = esstimating_fm_in_conv_circle(conv_radius(jj,1), sigma, lambda_in(idxx2), 'left');
                        
                        if ( isempty(idxx2) )
                            conv_region_angle = ang_thm_pt; %max(angle(tmp(idx_tmp)));
                            fprintf('conv_region_angle_1 = %11.4e\n', conv_region_angle);
                        else
                            if ( ((prev_no_ew_isempty >= 3) || ...
                                    (prev_no_ew_isempty >= 2 && mod(angle(shift(jj,1)),2*pi) > 100/180*pi) ...
                                    || (prev_no_ew_isempty >= 2 && length(idx3) >= 15)) ...
                                    && isempty(idx_neg_real_ew) )
                                
                                %weighting_choosing = 'left';
                                ang_cong_ew        = 0;
                                
                            else
                                
                                ang_cong_ew       = max(angle(lambda_in(idxx2)));
                                
                            end
                            
%                             idx5              = find(abs(lambda_in(idxx2)) < rad_old, 1);
%                             lambda_in(idxx2)
%                             if ( jj == 1 && isempty(idx5) )
%                                 ang_cong_ew   = 0;
%                             end
                            %ang_thm_pt        = max(angle(tmp(idx_tmp)));
                            %ang_cong_ew       = 0; %max(angle(lambda_in(idxx2)));
                            conv_region_angle = max(ang_thm_pt, ang_cong_ew);
                            
                            fprintf('conv_region_angle_2 = %11.4e with ang_thm_pt = %11.4e and ang_cong_ew = %11.4e \n', conv_region_angle, ang_thm_pt, ang_cong_ew);
                        end
                        
                        if ( isempty(idxx2) )
                            if ( jj == 1 )
                                theta_min = angle(shift(jj,1)) * 3;
                                fprintf('(1) theta_min = %11.4e \n', theta_min);
                                new_sigma = rad * exp(1i * theta_min);
                            elseif ( ~isempty(idx_neg_real_ew) )
                                [ new_sigma, shrink, arrive_right_tag, arrive_left_tag ] = check_shrinking_cond( ...
                                     ew_inside_c_new(no_ew_in_c_org+1:no_ew_in_c-mm2,1), rad, min_rad_neg_real_ew, ...
                                     sigma, conv_radius(jj,1), arrive_right_tag, arrive_left_tag, increment, lambda_in(idxx) );
                            else
                                right_angle = angle(shift(jj,1));
                                theta_min   = min(pi-2.0e-2, (right_angle+0.8*abs(right_angle - angle(shift(jj-1,1)))));
                                fprintf('(2) theta_min = %11.4e \n', theta_min);
                                new_sigma   = rad * exp(1i * theta_min); 
                            end
                        elseif ( ~isempty(idx_neg_real_ew) ) 
                            [ new_sigma, shrink, arrive_right_tag, arrive_left_tag ] = check_shrinking_cond( ...
                                  ew_inside_c_new(no_ew_in_c_org+1:no_ew_in_c-mm2,1), rad, min_rad_neg_real_ew, ...
                                  sigma, conv_radius(jj,1), arrive_right_tag, arrive_left_tag, increment, lambda_in(idxx) );
                        else
                            new_left_ew     = lambda_in(idxx2);
                            %right_angle     = max(angle(lambda_in(idxx2)));
                            %theta_min       = min(pi, (right_angle+abs(right_angle - angle(shift(jj,1)))/2));
                            idx3            = find(abs(lambda_in) < rad_old); 
                            if ( (prev_no_ew_isempty >= 3) || (prev_no_ew_isempty >= 2 && mod(angle(shift(jj,1)),2*pi) > 130/180*pi) ...
                                    || (prev_no_ew_isempty >= 2 && length(idx3) >= 15) || strcmp(weighting_choosing,'left') )
                                diff_angle         = abs(conv_region_angle - angle(shift(jj,1)));
                                                                
                                angle_ew       = angle(lambda_in(idxx2));
                                iddx           = find( abs(angle_ew) < 1.0e-9 );
                                if ( ~isempty(iddx) )
                                    angle_ew(iddx) = 0;
                                end
                                angle_ew       = mod(angle_ew,2*pi);
                                    
                                if ( diff_angle > 1.0e-5 )
                                    if ( (max(angle_ew) - min(angle_ew)) / diff_angle < 0.15 && jj == 1 )
                                        delta_theta = diff_angle*0.1;
                                    else
                                        delta_theta = min(abs(pi - conv_region_angle)*0.3, diff_angle*0.3);
                                    end
                                else
                                    delta_theta    = min(abs(pi - conv_region_angle)*0.3, abs(conv_region_angle - min(angle_ew))*0.8);
                                end
                                delta_theta        = min(10/180*pi, delta_theta);
                                weighting_choosing = 'left';
                                TT1         = [prev_no_ew_isempty delta_theta abs(pi - conv_region_angle)*0.3 abs(conv_region_angle - angle(shift(jj,1)))*0.3 10/180*pi]
                            else
                                delta_theta = min(abs(pi - conv_region_angle)/2, abs(conv_region_angle - angle(shift(jj,1)))/2); 
                                delta_theta = min(15/180*pi, delta_theta); 
                            end
                            theta_min       = min(pi-2.0e-2, conv_region_angle+delta_theta);
                            fprintf('(3) theta_min = %11.4e \n', theta_min);
                            new_sigma       = rad * exp(1i * theta_min); 
                            
                            idx_tmp = find(abs(ref_ew_inside_c + 1./ref_ew_inside_c - mu_0) <= conv_radius(jj,1), 1);
                            if ( ~isempty(idx_tmp) && abs(mu_0 + rad_old + 1/rad_old) <= conv_radius(jj,1) ...
                                    && abs(mu_0 + rad + 1/rad) <= conv_radius(jj,1) )
                                shrink     = 'yes';
                            end
                             
                        end
                        
                        shift(jj+1,1)    = new_sigma;
                        shift_1(jj+1, 1) = (shift(jj+1, 1)+1)./(shift(jj+1, 1)-1); 
                        
                    else
                          
                        new_left_ew  = lambda_in(idxx);
                            
                        if ( ~isempty(idx_neg_real_ew) )
                            increment = 'yes';
                            
                            [new_sigma, shrink, arrive_right_tag, arrive_left_tag] = check_shrinking_cond( ...
                                     ew_inside_c_new(no_ew_in_c_org+1:no_ew_in_c-mm2,1), rad, min_rad_neg_real_ew, ...
                                     sigma, conv_radius(jj,1), arrive_right_tag, arrive_left_tag, increment, lambda_in(idxx) );
                        else
                            
                            [new_sigma, new_left_ew, conv_region_angle, max_wanted_ew, m2] = choose_new_shift( sigma, lambda_in, ... %ew_new(idx_ew), ...
                                 ew_inside_c_new(1:no_ew_in_c,1), rad, rad_old, conv_radius(jj,1), max_wanted_ew_org, rat_max_angle, prev_no_ew_isempty );
                                                     
                             if ( isempty(new_left_ew) )
                                new_left_ew = lambda_in(idxx);
                             end
                        
                             if (angle(new_sigma) > 3.1067) % which is equal to 178/180*pi
                                 if (abs(min_rad_neg_real_ew + 1/min_rad_neg_real_ew - ...
                                            (sigma + 1 / sigma)) < conv_radius(jj,1))
                                     gamma_real = -rad * (1 + 1.0e-3); 
                                     gamma_imag = sqrt(gamma_real^2 - rad^2);
                                     new_sigma  = gamma_real + 1i * gamma_imag;
                                     %shrink     = 'yes';
                                     fprintf('shrink 3 with min_rad_neg_real_ew = %11.4e\n', min_rad_neg_real_ew);
                                 else
                                     idx_tmp   = find(angle(lambda_in(idxx)) >= pi - 3.0e-2 ...
                                             & angle(lambda_in(idxx)) < pi + 1.0e-3) ;
                                     if ( ~isempty(idx_tmp) )
                                         max_val   = max(abs(lambda_in((idxx(idx_tmp))))); 
                                         new_sigma = (min_rad_neg_real_ew - 2*max_val) / 3 + 1i * 1.0e-3 * rad;  
                                     end
                                     fprintf('shrink_next again 4 with new_sigma = %11.4e+1i(%11.4e)\n', real(new_sigma), imag(new_sigma));
                                 end
                             end
                        end
                        
                        shift(jj+1,1)    = new_sigma;
                        shift_1(jj+1, 1) = (shift(jj+1, 1)+1)/(shift(jj+1, 1)-1);
                        
                    end
                    disp(new_sigma)
                    
%                     if ( strcmp(arrive_right_tag,'yes') && abs(mod(angle(sigma),2*pi)) <= 179/180*pi )
%                         [ amplitude ] = esstimating_neg_am_in_conv_circle(conv_radius(jj,1), sigma);
%                         tmp_sigma     = amplitude * exp(1i * abs(mod(angle(sigma),2*pi)));
%                         disp(tmp_sigma)
%                     end
                                
                %
                % Detect shift value is outside the unit circle or not?
                % If yes, move it on the unit circle and stop it at next iteration
                % by setting flag_next_stop = 'yes'.
                %
            
                %
                % Find the candiate of the radius for the next loop
                %
                index2      = find(imag(lambda_in) >= -1.0e-11); %find(imag(ew_new)>=0);
                
                [ amplitude ] = esstimating_am_in_conv_circle(conv_radius(jj,1), sigma);
                
                if ( isempty(index2) || skill_flag == 1 )
                    ew_inner(jj, 1) = 0; %rad_old; 
                elseif ( skill_flag == 0 ) 
                    tmp    = lambda_in(index2); %ew_new(index);
                    index1 = find(abs(tmp) < rad);
                    if ( isempty(index1) )
                        ew_inner(jj, 1)     = 0; %rad_old; 
                        prediction_pt(jj,1) = amplitude * exp(1i * angle(sigma));
%                         if ( ~isempty(idx_tmp) )
%                             [~, idx_min]        = min( abs(tmp(idx_tmp)) );
%                             prediction_pt(jj,1) = tmp(idx_tmp(idx_min));
%                         end
                    else
                        [~,ij]          = min(abs(tmp(index1)));
                        ew_inner(jj, 1) = tmp(index1(ij)); %min(abs(tmp(index1))); 
                    end 
                    
                    if ( isempty(idx_ew) )
                        idx_rad_check = find(abs(ew_inside_c_new(no_ew_in_c-no_ref_ev+1:no_ew_in_c,1) - ew_inner(jj, 1)) < 1.0e-10, 1);
                        if ( ~isempty(idx_rad_check) )
                            ew_inner(jj, 1) = 0;
                        end 
                    end
                end  
                
% %  
% %             
% %                 rem = ['(' num2str(outer_it) ',' num2str(percentage) ')'];
% %                 figure(2)
% %                 if ( real(sigma) <= 0 )
% %                     text(real(sigma)*0.98,imag(sigma),rem);
% %                 else
% %                     text(real(sigma)*1.02,imag(sigma),rem);
% %                 end
            
 
                if ( strcmp(plot_flag,'yes') )
                    if ( skill_flag == 0 )
                        plot_ew_stepBystep( jj, ew_new, sigma, NO_ew_T, ew_new_1, shift, ew_skew_Cayley );
                
                        currFrame = getframe;
                        writeVideo(vidObj,currFrame);
                    else
                        figure(2) 
                        plot(real(sigma),imag(sigma),'kd','LineWidth',2); 
                        drawnow
                        hold on
                    
                        currFrame = getframe;
                        writeVideo(vidObj,currFrame);
        
                        max_wanted_ew = 10;
                    end
                end
                
                sigma = shift(jj+1,1);
                jj    = jj + 1;
    
                if ( strcmp(arrive_left_tag,'yes') || strcmp(arrive_right_tag,'yes') || ...
                        (angle(shift(jj,1)) >= 3.1154 && rad <= 0.05 ) ) % which is equal to 178.5/180*pi 
                    if ( rad > 0.05 )
                        max_no_re_estimate    = 0;
                        max_no_skew_J_Lan     = 1;
                    else
                        max_no_re_estimate    = -1;
                        max_no_skew_J_Lan     = 0;  
                    end
                    
                    deflation_enforce = 'no';
                    
                end
            else % GTSHIRA does not convergent within the maximal iteration
                sigma       = rad * exp(1i*min(179.5/180*pi,angle(sigma)+15/180*pi));
                shift(jj,1) = sigma;
            end
    
    end
end

if ( strcmp(shrink, 'yes') )
                            
    %tmp_vec               = [ref_ew_inside_c; ew_inside_c_new(no_ew_in_c-no_ref_ev+1:no_ew_in_c,1)]; % no_ew_in_c_org+1
    tmp_vec               = [ref_ew_inside_c; ew_inside_c_new(no_ew_in_c_org+1:no_ew_in_c,1)];
    idx_neg_real_ew       = find(angle(tmp_vec) >= pi - 3.0e-2 ...
                & angle(tmp_vec) < pi + 1.0e-3 & abs(tmp_vec) <= rad_old); % >= rad); %(rad + rad_old)/2);
    min_rad_neg_real_ew   = min(abs(tmp_vec(idx_neg_real_ew))); 
            
    idx_pos_real_ew       = find(angle(tmp_vec) >= - 1.0e-3 ...
                & angle(tmp_vec) < 3.0e-2 & abs(tmp_vec) <= rad_old); % >= rad);  %(rad + rad_old)/2);
    min_rad_pos_real_ew   = min(abs(tmp_vec(idx_pos_real_ew)));
    
    no_pos_real_ew        = length(idx_pos_real_ew);
    no_neg_real_ew        = length(idx_neg_real_ew);
            
    fprintf('number of negative real eigenvalues = %2.0f \n', no_neg_real_ew)
    fprintf('number of positive real eigenvalues = %2.0f \n', no_pos_real_ew)
               
%             idx_2 = abs(ew_inside_c_new(1:no_ew_in_c,1)) <= rad_old;
            
    diff_rad = abs(rad_old - rad);
    rad_old  = rad;
    
    [ rad, max_wanted_ew, radius_case ] = determine_radius_circle( ew_inner, rad_old, no_pos_real_ew, ...
        no_neg_real_ew, min_rad_pos_real_ew, min_rad_neg_real_ew, ...
        ew_inside_c_new(no_ew_in_c-no_ref_ev+1:no_ew_in_c,1), diff_rad, ...
        ew_inside_c_new(1:no_ew_in_c,1), prediction_pt(1:jj-1,1) );
    
    fprintf('Radius %8.5f shrink to %8.5f \n\n', rad_old, rad);  
    ref_ev_outside_c      = next_ref_ev_outside_c(:,1:no_ref_ev); 
    ref_ew_inside_c       = ew_inside_c_new(no_ew_in_c-no_ref_ev+1:no_ew_in_c,1);
    ref_ev_inside_c       = ev_inside_c_new(:,no_ew_in_c-no_ref_ev+1:no_ew_in_c);
       
    ew_inner              = ew_inner(1:jj,1);
end

if ( strcmp(checking_choice, 'yes') )
    if ( min(find(mod(angle(ref_ew_inside_c), 2*pi) >= 165/180*pi)) <= 178/180*pi )
        weighting_choosing = 'right';
    end
end

ew_inside_c_new  = ew_inside_c_new(1:no_ew_in_c,1);
ev_inside_c_new  = ev_inside_c_new(:,1:no_ew_in_c);
no_deuplicate_ew = no_deuplicate_ew(1:jj,:);

flag_conv        = 0;

end

