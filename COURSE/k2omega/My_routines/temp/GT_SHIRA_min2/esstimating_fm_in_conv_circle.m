function [ frequency ] = esstimating_fm_in_conv_circle(conv_rad, sigma, lambda, left_right_angle)
%
% Compute the solution x of the following nonlinear equation
%
%  ((x+1/x)*cos(theta)-(beta+1/beta)*cos(phi))^2 +
%  ((x-1/x)*sin(theta)-(beta-1/beta)*sin(phi))^2 - conv_rad^2 = 0
%
% where x = beta
%

    tol_angle = 1.0e-7;

    beta   = abs(sigma);
    phi    = angle(sigma);
    if ( abs(phi) <= tol_angle)
        phi = 0;
    end
    phi    = mod(phi, 2*pi);
    
    delta1 = ( beta + 1 / beta ) * cos(phi);
    delta2 = ( beta - 1 / beta ) * sin(phi);
    
    tmp_ang1  = angle(lambda);
    idxx      = find( abs(tmp_ang1) <= tol_angle );
    if ( ~isempty(idxx) )
        tmp_ang1(idxx) = 0;
    end
    
    %angle_sigma = mod(angle(sigma),2*pi);
    
    switch left_right_angle
        case 'left'
            if ( isempty(tmp_ang1) )
                initial_ang1 = phi*1.05;
            else
                initial_ang  = min(pi, max(mod(tmp_ang1,2*pi)));
                initial_ang1 = max(initial_ang, phi*1.05);
            end
            if ( initial_ang1 <= 5/180*pi )
                initial_ang1 = initial_ang1 * 2;
            end
            max_ratio    = 1.05;
            ratio_increm = 1.005;
            %fprintf(' left case with initial_ang1 = %11.4e \n',initial_ang1);
    
        case 'right'
            if ( isempty(tmp_ang1) )
                initial_ang1 = phi*0.95;
            else
                initial_ang  = min(mod(tmp_ang1,2*pi));
                initial_ang1 = min(initial_ang, phi*0.95);
            end
            max_ratio    = 0.95;
            ratio_increm = 0.995;
            %fprintf(' right case with initial_ang1 = %11.4e \n',initial_ang1);
            
    end

    theta0       = initial_ang1; %phi * 1.2;
    
    flag = 1;
    kk   = 1;
    while ( kk <= 3 && flag == 1 )
        [frequency, ~, ~, flag] = newton_method_funhandle(@(theta)fun_f_theta(beta, theta, delta1, delta2, conv_rad), ...
            @(theta)fun_df_theta(beta, theta, delta1, delta2), theta0, 1.0e-12);
        
        frequency = mod(frequency,2*pi);
        
        switch left_right_angle
            case 'left'
                if ( frequency <= phi || (frequency > pi+1.0e-2 && phi < 100/180*pi))
                    flag = 1; 
                    if ( frequency > pi+1.0e-2 )
                        %theta0    = max(min(initial_ang1*max_ratio^kk, theta0 * ratio_increm), (2*pi-frequency)* ratio_increm^kk);
                        theta0    = min(initial_ang1*max_ratio^kk, theta0 * ratio_increm) + (2*pi-frequency)* ratio_increm^kk; 
                    elseif ( frequency <= phi )
                        theta0    = max(min(initial_ang1*max_ratio, theta0 * ratio_increm), 2*phi-frequency);
                    end
                end
                
            case 'right'
                if ( frequency >= phi )
                    flag = 1;
                end
                theta0    = min(initial_ang1*max_ratio, theta0 * ratio_increm);
                
        end
        val       = fun_f_theta( beta, frequency, delta1, delta2, conv_rad );
        %fprintf('frequency(%1.0f) = %11.4e, angle_sigma = %11.4e, rsdl = %11.4e, initial = %11.4e  \n',kk, frequency, phi, val, theta0 );
        kk     = kk + 1;
    end
    
    switch left_right_angle
        case 'left'
            if ( flag == 1 || frequency < phi )% || frequency > initial_ang1*2.5) 
                frequency = initial_ang;
            end
            frequency = max(frequency, phi);
    
        case 'right'
            if ( flag == 1 || frequency > phi || mod(frequency, 2*pi) > pi || frequency < 0 )
                frequency = initial_ang; 
            end
            
    end 
    
end

function [ f_theta ] = fun_f_theta( beta, theta, delta1, delta2, rad )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

t1 = ( beta + 1/beta ) * cos(theta) - delta1;
t2 = ( beta - 1/beta ) * sin(theta) - delta2;

f_theta = t1^2 + t2^2 - rad^2;

end

function [ fd_theta ] = fun_df_theta(beta, theta, delta1, delta2)

t3 = beta + 1 / beta;
t4 = beta - 1 / beta;
t1 = t3 * cos(theta) - delta1;
t2 = t4 * sin(theta) - delta2;

fd_theta =  2 * t2 * t4 * cos(theta) - 2 * t1 * t3 * sin(theta);

end