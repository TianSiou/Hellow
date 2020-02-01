function [ new_sigma, shrink, arrive_right_tag, arrive_left_tag ] = check_shrinking_cond( ...
    ew_inside_c_prev, rad, min_rad_neg_real_ew, sigma, conv_radius, arrive_right_tag, ...
    arrive_left_tag, increment, current_lambda )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if ( strcmp(arrive_right_tag,'yes') && abs(mod(angle(sigma),2*pi)) > 178.5/180*pi )
    [ amplitude ] = esstimating_neg_am_in_conv_circle(conv_radius, sigma);
    tmp_sigma     = amplitude * exp(1i * abs(mod(angle(sigma),2*pi)));
    fprintf('Theoretical intersction point = %11.4e + (%11.4e)*1i \n', real(tmp_sigma), imag(tmp_sigma)); 
else
    tmp_sigma     = 0;
end
                    
shrink      = 'no';
ew_inside_c = [ew_inside_c_prev; current_lambda];

idx    = find(abs(imag(ew_inside_c)) <= 1.0e-8 & real(ew_inside_c) <= -1.0e-7);

if ( angle(sigma) >= 165/180*pi )
    if ( strcmp(arrive_right_tag,'no') )
        % The target eigenvalues in the right region are not computed in the
        % previous step. Check it compute in this step or not.
        if ( abs(-rad - 1/rad - sigma - 1/sigma) <= conv_radius || ...
                ~isempty(find(real(ew_inside_c(idx)) >= -rad, 1)) );
            arrive_right_tag = 'yes';
        end
    end

    if ( strcmp(arrive_left_tag,'no') )
        % The target eigenvalues in the left region are not computed in the
        % previous step. Check it compute in this step or not. 
        
        if ( abs(min_rad_neg_real_ew + 1/min_rad_neg_real_ew - sigma - 1/sigma) <= conv_radius || ...
                ~isempty(find(real(ew_inside_c(idx)) <= min_rad_neg_real_ew, 1)) || abs(min_rad_neg_real_ew) < rad );
            arrive_left_tag = 'yes'; 
            
        end
    end
elseif ( ~isempty(idx) )
    if ( strcmp(arrive_right_tag,'no') )
        % The target eigenvalues in the right region are not computed in the
        % previous step. Check it compute in this step or not.
        if ( abs(-rad - 1/rad - sigma - 1/sigma) <= conv_radius || ...
                ~isempty(find(real(ew_inside_c(idx)) >= -rad, 1)) );
            arrive_right_tag = 'yes';
        end
    end

    if ( strcmp(arrive_left_tag,'no') )
        % The target eigenvalues in the left region are not computed in the
        % previous step. Check it compute in this step or not.  
        
        if ( abs(min_rad_neg_real_ew + 1/min_rad_neg_real_ew - sigma - 1/sigma) <= conv_radius || ...
                ~isempty(find(real(ew_inside_c(idx)) <= min_rad_neg_real_ew, 1)) || abs(min_rad_neg_real_ew) < rad);
            arrive_left_tag = 'yes';
        end
    end
end

% idx_lambda = find(abs(imag(lambda)) <= 1.0e-8 & real(lambda) <= -rad);
% 
if ( strcmp(increment,'yes') )
    
     %weight  = min(1,abs(max(real(ew_inside_c(idx))) - min(ew_inside_c(idx))) / abs(min_rad_neg_real_ew + rad) * 1.7);
     if ( strcmp(arrive_right_tag,'yes') )
         % all the target eigenvalues in the right region are computed
         % goto the left region
         % 
         if ( strcmp(arrive_left_tag,'yes') )
             % all the target eigenvalues in the left region are also computed
             %
             shrink    = 'yes';
             new_sigma = sigma;
             fprintf('shrink I in check_shrinking_cond \n');
         else
             % Move the shift value from min(real(current_lambda))) to the left region
             % 
             
             idx_tmp       = find(abs(imag(ew_inside_c_prev)) <= 1.0e-8 & real(ew_inside_c_prev) <= -1.0e-7);
             [~,idx_tmp1]  = max(real(ew_inside_c_prev(idx_tmp)));
             right_ew_prev = ew_inside_c_prev(idx_tmp(idx_tmp1));
             idx_tmp2      = imag(ew_inside_c) > -1.0e-8 & imag(ew_inside_c) < rad*10 & real(ew_inside_c) <= -1.0e-7;
             angle_tmp     = min(179.5 / 180 * pi, max(177 / 180 * pi, (pi + 2 * max(angle(ew_inside_c(idx_tmp2)))) / 3));
             
             if ( isempty(idx_tmp1) )
                 min_left_ew = min(real(ew_inside_c(idx)));
                 delta_sigma = abs(min(-rad, real(sigma)) - min_left_ew);
                 if ( delta_sigma ~= 0 )
                     ratio_sigma = abs(min_rad_neg_real_ew - min_left_ew) / delta_sigma;
                 else
                     ratio_sigma = 10^7;
                 end
                 if ( ratio_sigma >= 0.95 )
                     delta_sigma = 0.4 * delta_sigma;
                 else
                     delta_sigma = 0.65 * delta_sigma;
                 end 
                 new_sigma   = min_left_ew - delta_sigma;
                 if ( new_sigma > min_rad_neg_real_ew )
                     if (abs(new_sigma-min_rad_neg_real_ew) < abs(new_sigma-min_left_ew))
                         new_sigma = new_sigma + (abs(new_sigma-min_left_ew) - abs(new_sigma-min_rad_neg_real_ew));
                     end
                 else
                     new_sigma = (2*min_left_ew + min_rad_neg_real_ew) / 3;
                 end 
                 new_sigma = abs(new_sigma) * exp(1i * angle_tmp); 
                 fprintf('(a0) shrink_next 1 with sigma = %11.4e+1i(%11.4e) in checking_shrinking_cond \n', real(new_sigma), imag(new_sigma)); 
             elseif ( abs(right_ew_prev+1/right_ew_prev-(sigma+1/sigma)) < conv_radius || real(right_ew_prev) > real(sigma))
                 %
                 % revision at 09/14
                 %
                 idx_lambda  = real(current_lambda) < 0 & imag(current_lambda) > -1.0e-11 & imag(current_lambda) < 1.0e-7;
                 min_left_ew = min(real(current_lambda(idx_lambda)));
                 
                 %min_left_ew = min(real(ew_inside_c(idx))); original version
                 idx_neg_ew = find(imag(current_lambda) > -1.0e-11 & imag(current_lambda) < 1.0e-7 & real(current_lambda) <= -rad );
                 if ( ~isempty(idx_neg_ew) )
                     dist_ew    = abs(min_left_ew - max(real(current_lambda(idx_neg_ew))));
                     density_ew = length(idx_neg_ew) / dist_ew;
                 end
                 
                 delta_sigma = abs(min(-rad, real(sigma)) - min_left_ew);
                 if ( delta_sigma ~= 0 )
                     ratio_sigma = abs(min_rad_neg_real_ew - min_left_ew) / delta_sigma;
                 else
                     ratio_sigma = 10^7;
                 end
                 
                 %[ rad real(sigma) min_left_ew min_rad_neg_real_ew delta_sigma ratio_sigma density_ew ]
                 
                 %if ( ratio_sigma >= 2 )
                 %    delta_sigma = 0.35 * delta_sigma;
                 %else
                 if ( ratio_sigma >= 1.1 && ~isempty(idx_neg_ew) )
                     delta_sigma = min(0.4, 100/sqrt(density_ew)) * delta_sigma;
                 elseif ( ratio_sigma >= 0.95 )
                     delta_sigma = 0.4 * delta_sigma;
                 else
                     delta_sigma = 0.65 * delta_sigma;
                 end
                 %delta_sigma = min(delta_sigma, 0.95 * new_leng_interval / 2);
                 new_sigma   = min_left_ew - delta_sigma;
                 if ( new_sigma > min_rad_neg_real_ew )
                     if (abs(new_sigma-min_rad_neg_real_ew) < abs(new_sigma-min_left_ew))
                         new_sigma = new_sigma + (abs(new_sigma-min_left_ew) - abs(new_sigma-min_rad_neg_real_ew));
                         fprintf('(i)\n');
                     end
                     fprintf('(ii)\n');
                 else
                     new_sigma = (2*min_left_ew + min_rad_neg_real_ew) / 3; 
                     fprintf('(iii)\n');
                 end
%              new_sigma   = max(min_rad_neg_real_ew, min(real(ew_inside_c(idx))) - 0.8*abs(min(-rad,real(sigma)) ...
%                  - min(real(ew_inside_c(idx))))); 
                 new_sigma = abs(new_sigma) * exp(1i * angle_tmp);
%              new_sigma = max(min_rad_neg_real_ew, min(real(ew_inside_c(idx))) - 0.8*abs(min(-rad,real(sigma)) ...
%                  - min(real(ew_inside_c(idx))))) + 1i * 1.0e-3 * rad;
                 fprintf('(a) shrink_next 1 with sigma = %11.4e+1i(%11.4e) in check_shrinking_cond \n', real(new_sigma), imag(new_sigma)); 
             else 
                 idx_tmp       = find(abs(imag(current_lambda)) <= 1.0e-8 & real(current_lambda) <= -1.0e-7);
                 if ( ~isempty(idx_tmp) )
                     [min_val,idx_tmp1] = min(real(current_lambda(idx_tmp)));
                     max_val            = max(real(current_lambda(idx_tmp)));
                     if ( length(idx_tmp) == 1 )
                         ratio_leng     = abs(max_val - rad) / abs(right_ew_prev - min_val);
                     else
                         %density            = length(idx_tmp) / abs(max_val - min_val);
                         ratio_leng     = abs(max_val - min_val) / abs(right_ew_prev - min_val);
                     end
                     left_ew_cur        = current_lambda(idx_tmp(idx_tmp1));
                     %ratio_leng
                     if ( ratio_leng >= 0.8 )
                         new_sigma     = (left_ew_cur + right_ew_prev) / 2 + 1.0e-3*rad; %* exp(1i * 179.3 / 180 * pi);
                     else
                         tt            = ratio_leng / 2;
                         new_sigma     = (1 - tt) * left_ew_cur + tt * right_ew_prev + 1.0e-3*rad;
                     end
                     %[left_ew_cur right_ew_prev]
                 else
                     new_sigma     = right_ew_prev + 1.0e-3*rad;
                 end
                 
                 fprintf('(f) shrink_next 1 with sigma = %11.4e+1i(%11.4e) in checking_shrinking_cond \n', real(new_sigma), imag(new_sigma));
             end
             
         end
    
     else 
         if ( strcmp(arrive_left_tag,'yes') )
             % all the target eigenvalues in the left region are computed
             % Move the shift value from max(real(ew_inside_c(idx)))) to the right region
             %
             %new_sigma   = ( (1 - weight) * rad + weight * max(real(ew_inside_c(idx)))) + 1i * 1.0e-3 * rad;
             max_right_ew = max(real(ew_inside_c(idx))); 
             delta_sigma  = abs( max_right_ew - min(-rad,real(sigma)) );
             
             if ( delta_sigma ~= 0 )
                 ratio_sigma = abs(max_right_ew - min(-rad,real(sigma))) / delta_sigma;
             else
                 ratio_sigma = 10^7;
             end
             
             if ( ratio_sigma >= 0.95 )
                 delta_sigma = 0.4 * delta_sigma;
             else
                 delta_sigma = 0.65 * delta_sigma;
             end
             
             new_sigma   = min(-rad, max(real(ew_inside_c(idx))) + delta_sigma);
             new_sigma   = abs(new_sigma) * exp(1i * 177 / 180 * pi);
%              new_sigma   = min(-rad, max(real(ew_inside_c(idx))) + 0.65*abs(max(real(ew_inside_c(idx))) - min(-rad,real(sigma)))) + 1i * 1.0e-3 * rad;
             fprintf('(b) shrink_next 2 with sigma = %11.4e+1i(%11.4e) in checking_shrinking_cond \n', real(new_sigma), imag(new_sigma));
         else
             % Compute the shift value in the left region
             % 
             %new_sigma   = min(-rad, max(real(ew_inside_c(idx))) + 0.65*abs(max(real(ew_inside_c(idx))) - min(-rad,real(sigma)))) + 1i * 1.0e-3 * rad;
             new_sigma = rad * exp(1i * 178 / 180 * pi);
             fprintf('(e) shrink_next 3 with sigma = %11.4e+1i(%11.4e) in checking_shrinking_cond \n', real(new_sigma), imag(new_sigma));
         end
    
    end

else
    %
    % All the eigenvalues are duplicated. No new nonegative eigenvalues are added to the set idx
    % The rule of choosing the shift value needs to modify
    %
    if ( strcmp(arrive_right_tag,'yes') )
        
        if ( strcmp(arrive_left_tag,'yes') )
             % all the target eigenvalues in the left region are also computed
             %
             shrink    = 'yes';
             new_sigma = sigma;
             fprintf('shrink II in check_shrinking_cond \n');
        else
            % Move the shift value more near to the left point
            % 
            new_sigma = 1.05 * real(sigma) + 1i * 1.0e-3 * rad;
            %weight = min(1, abs(max(real(ew_inside_c(idx))) - min(ew_inside_c(idx))) / abs(min_rad_neg_real_ew + rad) / 1.3);
            %     new_sigma = ( weight * min_rad_neg_real_ew + (1 - weight) * min(real(ew_inside_c(idx)))) + 1i * 1.0e-3 * rad;
            fprintf('(c) shrink_next 3 with sigma = %11.4e+1i(%11.4e) in checking_shrinking_cond \n', real(new_sigma), imag(new_sigma));
        end
    
     else 
%          if ( strcmp(arrive_left_tag,'yes') )
             % all the target eigenvalues in the left region are computed
             % Move the shift value more near to the right point
             %
        %weight = min(1, 1.2 * abs(max(real(ew_inside_c(idx))) - min(ew_inside_c(idx))) / abs(min_rad_neg_real_ew + rad));
             new_sigma   = 0.95 * real(sigma) + 1i * 1.0e-3 * rad;
             %new_sigma   = ( (1 - weight) * rad + weight * max(real(ew_inside_c(idx)))) + 1i * 1.0e-3 * rad;
             fprintf('(d) shrink_next 4 with sigma = %11.4e+1i(%11.4e) in checking_shrinking_cond \n', real(new_sigma), imag(new_sigma));
%          else
%              % Compute the shift value in the left region
%              % 
%              new_sigma   = ( min_rad_neg_real_ew / 3 + 2 / 3 * min(real(ew_inside_c(idx)))) + 1i * 1.0e-3 * rad;
%              fprintf('shrink_next 3 with sigma = %11.4e+1i(%11.4e) \n', real(new_sigma), imag(new_sigma));
%          end
    
    end
    
end

if ( strcmp(arrive_right_tag,'yes') && abs(tmp_sigma) > abs(new_sigma) )
    new_sigma = tmp_sigma;
end

end

