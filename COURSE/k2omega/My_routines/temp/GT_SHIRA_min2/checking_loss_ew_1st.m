function [ doing_GTSHIRA_again, shift ] = checking_loss_ew_1st( lambda, prev_ew, ref_ew_inside_c, ...
    sigma, conv_radius, rad, pre_rad, ew_no_dupli )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

tol_angle        = 1.0e-10;

idx_angle        = imag(lambda) >= -1.0e-10 & abs(lambda) <= 1;
right_ang_thm_pt = esstimating_fm_in_conv_circle(conv_radius, sigma, lambda(idx_angle), 'right');

idx_angle   = find(imag(lambda) >= -1.0e-10 & abs(lambda) <= 1); 

if ( ~isempty(idx_angle) )
    tmp_ang2       = angle(lambda(idx_angle));
    idxx           =  abs(tmp_ang2) <= 1.0e-10 ;
    tmp_ang2(idxx) = 0;
    tmp_ang2       = mod(tmp_ang2,2*pi); 
    right_angle    = min(right_ang_thm_pt, min(tmp_ang2)); 
else
    right_angle = right_ang_thm_pt;
end

ref_angle = angle(ref_ew_inside_c);
idxx      = find( abs(ref_angle) <= tol_angle );
if ( ~isempty(idxx) )
    ref_angle(idxx) = 0;
end
ref_angle = mod(ref_angle,2*pi);

idx1 = find( abs(prev_ew) <= pre_rad & imag(prev_ew) >= -tol_angle );
idx2 = find( abs(ref_ew_inside_c) <= pre_rad & ref_angle >= right_angle);

if ( isempty(idx1) && isempty(idx2) )
    doing_GTSHIRA_again = 'no';
    shift               = [];
else
    if ( isempty(idx1) ) 
        idx4 = find(abs(ref_ew_inside_c(idx2) + 1./ref_ew_inside_c(idx2) - (sigma + 1 / sigma)) <= conv_radius, 1);
        if ( ~isempty(idx4) )
            doing_GTSHIRA_again = 'no';
            shift               = [];
        else
            doing_GTSHIRA_again = 'yes';
            shift               = (min(abs(ref_ew_inside_c(idx2))) + max(abs(lambda))) / 2;
            if ( sigma < 0 )
                angle_shift         = (mod(angle(-sigma),2*pi) + right_angle) / 2;
            else 
                angle_shift         = max((3*mod(angle(sigma),2*pi) + 2*right_angle) / 5, 1/180*pi);
            end

            shift               = shift * exp(1i * angle_shift);
            fprintf('Loss top ew (a) with theta = %13.6e \n',angle_shift);
        end
    else
        idx3 = find(abs(prev_ew(idx1) + 1./prev_ew(idx1) - (sigma + 1 / sigma)) <= conv_radius, 1);
        if ( ~isempty(idx3) )
            if ( isempty(idx2) )
                doing_GTSHIRA_again = 'no';
                shift               = [];
            else
                idx4 = find(abs(ref_ew_inside_c(idx2) + 1./ref_ew_inside_c(idx2) - (sigma + 1 / sigma)) <= conv_radius, 1);
                if ( ~isempty(idx4) )
                    doing_GTSHIRA_again = 'no';
                    shift               = [];
                else
                    doing_GTSHIRA_again = 'yes';
                    shift               = (min(abs(ref_ew_inside_c(idx2))) + max(abs(lambda))) / 2;
                    angle_ew            = angle(prev_ew(idx1));
                    idx_ang             =  angle_ew >= 0 ; 
                    angle_ew            = angle_ew(idx_ang); 
                    angle_current_ew    = mod(angle(ew_no_dupli), 2*pi);
                    angle_shift         = (min(max(angle_current_ew),mod(angle(sigma),2*pi)) + 2*max(mod(angle_ew,2*pi))) / 3;
                    %angle_shift         = (mod(angle(sigma),2*pi) + 2*max(mod(angle_ew,2*pi))) / 3;
                    shift               = shift * exp(1i * angle_shift);
                    fprintf('Loss ew left (b) with theta = %13.6e \n',angle_shift); 
                end
            end
        else
            idx5 = find( abs(prev_ew(idx1)) >= rad, 1 );
            if ( isempty(idx5) )
                doing_GTSHIRA_again = 'no';
                shift               = [];
            else
                doing_GTSHIRA_again = 'yes';
                tmp_ang             = angle(prev_ew(idx1));
                idx_tmp_ang         = find(abs(tmp_ang) <= tol_angle);
                if ( ~isempty(idx_tmp_ang) )
                    tmp_ang(idx_tmp_ang) = 0;
                end
                tmp_ang             = mod(tmp_ang,2*pi); 
                angle_shift         = (right_angle + max(tmp_ang)) / 2;
                if ( ~isempty(idx2) )
                    shift           = (min(abs(ref_ew_inside_c(idx2))) + max(abs(lambda))) / 2;
                else
                    shift           = max(abs(lambda));
                end
                shift               = shift * exp(1i * angle_shift);
                fprintf('Loss left ew (c) with theta = %13.6e \n',angle_shift);
                fprintf('right_angle = %13.6e, max_ang_ew = %11.6e \n',right_angle, max(angle(prev_ew(idx1))));
            end
        end
    end


 
end

end

