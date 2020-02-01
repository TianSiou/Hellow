function [ amplitude ] = esstimating_neg_am_in_conv_circle(conv_rad, sigma)
%
% Compute the solution x of the following nonlinear equation
%
%  ((x+1/x)*cos(theta)-(beta+1/beta)*cos(phi))^2 +
%  ((x-1/x)*sin(theta)-(beta-1/beta)*sin(phi))^2 - conv_rad^2 = 0
%
% where theta = phi
%
    beta   = abs(sigma);
    theta  = angle(sigma);
    delta1 = ( beta + 1 / beta ) * cos(theta);
    delta2 = ( beta - 1 / beta ) * sin(theta);
    x0     = beta * 1.05;
    
    flag = 1;
    kk   = 1;
    while ( kk <= 3 && flag == 1 )
        [amplitude, ~, rsdl, flag] = newton_method_funhandle(@(x)fun_f_rad(x, theta, delta1, delta2, conv_rad), ...
            @(x)fun_df_rad(x, theta, delta1, delta2), x0, 1.0e-12);
        
%         if ( abs(amplitude) > 1 )
%             amplitude = 1/ amplitude;
%             rsdl      = fun_f_rad(amplitude, theta, delta1, delta2, conv_rad);
%         end
        %fprintf('%1.0f: x = %11.4e + (%11.4e)*1i with rsdl = %11.4e \n', kk, real(amplitude), imag(amplitude), rsdl);
        
        x0 = max(0.99, x0 * 1.05);
        kk = kk + 1;
    end
    
    if ( flag == 1 )
        amplitude = 0;
    end
    
end

function [ f_x ] = fun_f_rad( x, theta, delta1, delta2, rad )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

t1 = ( x + 1/x ) * cos(theta) - delta1;
t2 = ( x - 1/x ) * sin(theta) - delta2;

f_x = t1^2 + t2^2 - rad^2;

end

function [ fd_x ] = fun_df_rad(x, theta, delta1, delta2)

t1 = ( x + 1/x ) * cos(theta) - delta1;
t2 = ( x - 1/x ) * sin(theta) - delta2;

fd_x = 2 * t1 * ( 1 - x^(-2) ) * cos(theta) + 2 * t2 * ( 1 + x^(-2)) * sin(theta);

end