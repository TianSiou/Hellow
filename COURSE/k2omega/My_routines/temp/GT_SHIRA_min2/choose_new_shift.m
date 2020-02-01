function [ new_sigma, new_left_ew, right_angle, max_wanted_ew, mod_no ] = choose_new_shift( sigma, ...
    lambda, ew_inside_c_new, rad, prev_rad, conv_radius, max_wanted_ew_org, rat_max_angle, prev_no_ew_isempty )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Let
%      C = {x; | x - sigma - 1/sigma | <= conv_radius}
% be the convergent region of GTSHIRA with shift sigma.
%
% Find the angle 'right_angle' of x = rad * exp(1i*right_angle) satisfying
%               | x - sigma - 1/sigma | = conv_radius
% 
% density     = 0;
idx_angle   = find(imag(lambda) >= -1.0e-10 & abs(lambda) <= 1);

% mu_0        = sigma + 1 / sigma;

left_ang_thm_pt  = esstimating_fm_in_conv_circle(conv_radius, sigma, lambda(idx_angle), 'left');
right_ang_thm_pt = esstimating_fm_in_conv_circle(conv_radius, sigma, lambda(idx_angle), 'right');

amp         = abs(sigma);
theta       = acos(real(sigma / amp));

density_rat = 5; %75;
angle_rat   = 0.96;

if ( ~isempty(idx_angle) )
    tmp_ang2    = angle(lambda(idx_angle)); 
    min_tmp_ang = min(tmp_ang2);
    max_tmp_ang = max(tmp_ang2);
    right_angle = max(theta, max(angle_rat*left_ang_thm_pt, max_tmp_ang)); %max(theta, 0.96*max(left_ang_thm_pt, max_tmp_ang));
    min_angle   = min(right_ang_thm_pt, min_tmp_ang); %min(min(tmp_ang1), min(tmp_ang2));
    %fprintf('(1) right_angle = %11.4e and max(tmp_ang2) = %11.4e in choose_new_shift \n', right_angle, max(tmp_ang2))
else
    right_angle = angle_rat*left_ang_thm_pt; % 0.9 * left_ang_thm_pt;
    min_angle   = right_ang_thm_pt; %min(tmp_ang1);
    %fprintf('(2) right_angle = %11.4e in choose_new_shift  \n', right_angle)
end
   
%
% The angle of the new shift value 'new_sigma' must be greater than right_angle
%

if ( isempty(idx_angle) )  
    choose_angle  = min(pi-1.0e-2, theta*1.1);
    new_sigma     = amp * exp(1i*choose_angle);
    
    idx2          = find(real(ew_inside_c_new) <= real(new_sigma));
    [~,idx]       = sort(abs(new_sigma-ew_inside_c_new(idx2)));
    new_left_ew   = ew_inside_c_new(idx2(idx(1:10)));
    
    mod_no        = 8;
    max_wanted_ew = max_wanted_ew_org;
else
    
%     right_angle = max(angle(lambda(idx_angle)));
    
    ew_angle      = angle(ew_inside_c_new);
    idx_tmp       = find( abs(ew_angle) <= 1.0e-10 );
    if ( ~isempty(idx_tmp) )
        ew_angle(idx_tmp) = 0;
    end
    ew_angle = mod(ew_angle, 2*pi);
    
    idx_region_ew = find( ew_angle > min_angle & ew_angle < right_angle & abs(ew_inside_c_new) <= prev_rad);
    if ( length(lambda) == 1 )
        max_angle     = 12 / 180 * pi;
        mod_no        = 8;
        max_wanted_ew = max_wanted_ew_org;
    else 
        density     = length(idx_region_ew) / abs(right_angle - min_tmp_ang);% min_angle); 
        %[max_tmp_ang right_angle  min_tmp_ang]
        tt          = density_rat/density;
        if ( tt < 0.08 )
            tt = min(0.08, 2 * tt);
        end
        max_angle   = min(rat_max_angle/180*pi, tt); %min(17/180*pi, tt); %min_angle)); 
        %[rat_max_angle/180*pi  tt]
        fprintf('density of eigenvalues = %10.4e, max_angle = %10.4e \n', density, max_angle); 
        if ( density < 50 )
            mod_no        = 8;
            max_wanted_ew = max_wanted_ew_org;
        elseif ( density < 100 )
            mod_no        = 12;
            max_wanted_ew = min(25,max(20, max_wanted_ew_org + 7));
        elseif ( density < 400 )
            mod_no        = 18;
            max_wanted_ew = min(40,max(35, max_wanted_ew_org + 10));
        else
            mod_no        = 25;
            max_wanted_ew = min(50,max(45, max_wanted_ew_org + 15));
        end

    end
    choose_angle = max(max_tmp_ang*1.01, right_angle + max_angle);
    %[max_tmp_ang*1.01, right_angle + max_angle]
    
    if ( choose_angle > 165/180*pi )
        %fprintf('choose_angle = %11.4e \n', choose_angle); 
        if ( (prev_no_ew_isempty >= 3) || (prev_no_ew_isempty >= 2 && mod(angle(sigma),2*pi) > 130/180*pi ) )
            if ( length(lambda) >= 20 )
                choose_angle = min(pi-1.0e-2, (0.825*right_angle + 0.175*pi)); % min(pi-1.0e-2, (0.825*right_angle + 0.175*pi));
            else
                choose_angle = min(pi-1.0e-2, (4*right_angle + pi) / 5);
            end
        else
            choose_angle = min(pi-1.0e-2, (2*right_angle + pi) / 3);
        end
        %fprintf('new choose_angle = %11.4e \n', choose_angle);
    end
    new_sigma    = rad * exp(1i*choose_angle);
    diff_rad     = (prev_rad - rad) / 2;
    delta_theta  = diff_rad * 1.5;
    
    flag     = 0;
    kk       = 1;
    while ( kk <= 5 && flag == 0 && prev_rad < 1) % 5
        %
        % Check the new shift is located in the convergent eigenvalues or not.
        %
        idx = find( angle(ew_inside_c_new) > choose_angle & angle(ew_inside_c_new) <= min(pi,choose_angle+delta_theta) ...
            & abs(ew_inside_c_new) < abs(new_sigma));
        idx2 = find( angle(ew_inside_c_new) <= choose_angle & angle(ew_inside_c_new) > max(0,choose_angle-delta_theta) ...
            & abs(ew_inside_c_new) < abs(new_sigma)); 
        %fprintf('Left bottom = %2.0f, Left Top = %2.0f, diff_rad = %11.4e \n', length(idx), length(idx2), diff_rad);
        if ( length(idx) >= 1 && length(idx2) >= 1 )
            % the new shift is located in the convergent eigenvalues
            % Choose another shift
            
            %figure(2)
            %plot(real(new_sigma),imag(new_sigma),'k+','LineWidth',2);
            
            if ( angle(new_sigma) >= pi * 5 / 180 && angle(new_sigma) <= pi * 30 / 180 )
                choose_angle = choose_angle * 1.05;
                new_sigma    = amp * exp(1i*choose_angle); % * 0.98;
                kk           = kk + 1;
            elseif ( angle(new_sigma) <= pi * 0.95 )
                amp          = min(abs(ew_inside_c_new(idx))) * 0.98;
                choose_angle = choose_angle + (choose_angle - right_angle) / (2*kk+1);
                %choose_angle = choose_angle * 1.05;
                new_sigma    = amp * exp(1i*choose_angle);
                kk           = kk + 1;
            else
                flag = 1;
            end
        else
            flag         = 1;
        end
    end

    if ( kk > 1 )
        idx2        = find(real(ew_inside_c_new) <= real(new_sigma));
        [~,idx]     = sort(abs(new_sigma-ew_inside_c_new(idx2)));
        new_left_ew = ew_inside_c_new(idx2(idx(1:10)));
    else
        new_left_ew = [];
    end

end
end

