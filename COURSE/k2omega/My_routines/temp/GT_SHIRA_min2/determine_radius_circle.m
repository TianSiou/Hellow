function [ rad, max_wanted_ew, radius_case ] = determine_radius_circle( ew_inner, current_rad, no_pos_real_ew, ...
    no_neg_real_ew, min_rad_pos_real_ew, min_rad_neg_real_ew, all_ew_current_rad, diff_rad, ew_inside_c_new, ...
    prediction_pt, radius_moving )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%rad_ratio    = 1;

idx_rad      = find(abs(ew_inner) > 0);
idx_prev_ew  = find( abs(ew_inside_c_new) < current_rad);

if ( isempty(idx_rad) && isempty(idx_prev_ew) )
    radius_case   = 'normal';
    max_wanted_ew = 10;
    rad           = current_rad - 2 * diff_rad;
    if ( rad <= 0 )
        rad = current_rad * 0.95;
    end
else
    
    if ( isempty(idx_rad) )
        ew_inner = ew_inside_c_new(idx_prev_ew);
        idx_rad  = find(abs(ew_inner) > 0);
    end
    
    %ew_inner(idx_rad,1)
    
    [~,sort_idx] = sort(abs(ew_inner(idx_rad,1))); 

    ew_sort      = ew_inner(idx_rad(sort_idx),1);
    target_ew_1  = ew_sort(1);

    if ( length(sort_idx) > 1 )
    
        target_ew_2      = ew_sort(2);
        ratio_radius     = abs(ew_sort(2:end,1))./abs(ew_sort(1:end-1,1));
    
        idx_ratio_radius = find( ratio_radius > 1.15);
        type_cond_1      = 1;
        if ( ~isempty(idx_ratio_radius) )
            diff_angle = abs(mod(angle(ew_sort(idx_ratio_radius(1:end))), 2*pi) - ...
                         mod(angle(ew_sort(idx_ratio_radius(1:end)+1)), 2*pi));
            idx_diff_angle = find( diff_angle > 6/180*pi );
            if ( ~isempty(idx_diff_angle) )
                type_cond_1 = 0;
                target_ew_1 = ew_sort(idx_ratio_radius(idx_diff_angle(end)));
                target_ew_2 = ew_sort(idx_ratio_radius(idx_diff_angle(end))+1);
            end
        end
    
        idx_ratio_radius = find( ratio_radius > 1.114, 1);
        type_cond_2      = 1;
        anggle_1st_ew    = angle(ew_inner(idx_rad(1),1));
        if ( ~isempty(idx_ratio_radius) && anggle_1st_ew >= 160/180*pi )
            type_cond_2 = 0;
        end
        test_ratio    = abs(target_ew_2) / abs(target_ew_1);
    
        if ( type_cond_1 == 0 || type_cond_2 == 0 )
        %if ( abs(abs(ew_inner(idx_rad(sort_idx(2)),1)) - current_rad)/current_rad > 2.0e-2 )
            choosing_type = 1;
            radius_case   = 'weighting';
        else
            radius_case   = 'normal';
            choosing_type = 2;
        end
        
        no_ew = 5;
        if ( length(sort_idx) > no_ew )
            [~,sort_nidx] = sort(-abs(ew_inner(idx_rad,1))); 

            ew_sort2      = ew_inner(idx_rad(sort_nidx),1);
            abs_ew        = abs(ew_sort2);
            idx_target_ew = [ find( abs(abs_ew(1:end-1) - abs_ew(2:end)) > 1.0e-10*abs_ew(end) ); length(abs_ew) ];
            angle_all_ew   = mod(angle(ew_sort2(idx_target_ew,1)), 2*pi);
            [~,sort_angle] = sort(angle_all_ew);
            angle_all_ew   = angle_all_ew(sort_angle);
            if ( length(angle_all_ew) >= no_ew )
                if ( angle_all_ew(no_ew,1) - angle_all_ew(1,1) < 15/180*pi  )
                    choosing_type = 3;
                end
            end
%         elseif ( length(sort_idx) >= 3 )
%             ang_ew_sort = angle(ew_sort);
%             idx_ang     = find( abs(ang_ew_sort) <= 1.0e-7 );
%             if ( ~isempty(idx_ang) )
%                 ang_ew_sort(idx_ang) = 0;
%             end
%             ang_ew_sort        = mod(ang_ew_sort,2*pi);
%             prev_rad           = current_rad + diff_rad;
%             idx_ew_lt_prev_rad =  abs(ew_inside_c_new) <= prev_rad;
%             ew_lt_prev_rad     = ew_inside_c_new(idx_ew_lt_prev_rad);
        end
        
    else
        radius_case   = 'normal';
        choosing_type = 2;
        %rad_ratio     = 0.96;
    end

    switch choosing_type
        case 1
            if ( test_ratio >= 1.25 )
                %idx1 = find(abs(ew_inner(idx_rad,1) - ew_inner(idx_rad(sort_idx(1)),1)) <= 1.0e-12);
                idx1 = find(abs(ew_sort - target_ew_1) <= 1.0e-12);
                tt   = mod(angle(target_ew_1),2*pi);
                if ( length(idx1) == 1 )
                    left_ang1  = min(pi+0.001, tt + 15/180*pi);
                    right_ang1 = max(0, tt - 15/180*pi);
                else 
                    left_ang1  = min(pi+0.001, tt);
                    right_ang1 = max(0,tt - 30/180*pi);
                end
            
                ang_ew = mod(angle(all_ew_current_rad), 2*pi);
                idx1   = find( ang_ew >= right_ang1 & ang_ew <= left_ang1);
                if ( isempty(idx1) )
                    idx_all_ew =  abs(ew_inside_c_new) <= current_rad ;
                    ang_ew     = mod(angle(ew_inside_c_new(idx_all_ew)), 2*pi);
                    idx1       = find( ang_ew >= right_ang1 & ang_ew <= left_ang1);
                end
            
                idx2 = find(abs(ew_sort - target_ew_2) <= 1.0e-12);
                tt   = mod(angle(target_ew_2),2*pi);
                if ( length(idx2) == 1 )
                    left_ang2  = min(pi+0.001, tt + 15/180*pi);
                    right_ang2 = max(0, tt - 15/180*pi);
                else
                    left_ang2  = min(pi+0.001, tt);
                    right_ang2 = max(0, tt - 30/180*pi);
                end
            
                idx2   = find( ang_ew >= right_ang2 & ang_ew <= left_ang2);
             
                if ( length(idx2)/length(idx1) >= 7 )
                    weighting = max(0.08, length(idx1)/length(idx2)) * 3.5;
                    fprintf('(a) weighting = %11.4e and ratio = %11.4e in determine_radius_circle \n',weighting, length(idx1)/length(idx2));
                else
                    weighting = 0.5;
                    fprintf('(b) weighting = %11.4e in determine_radius_circle \n',weighting);
                end
            else
                weighting     = min(2 / 3 / abs(test_ratio - 0.114), 0.99);
                fprintf('(c) weighting = %11.4e in determine_radius_circle \n',weighting);
            end 
        
            rad           = (1-weighting) * abs(target_ew_2) + weighting * abs(target_ew_1); 
            max_wanted_ew = max(15, min(32, max(no_neg_real_ew, no_pos_real_ew)));
        
        case 2
        
            if ( length(sort_idx) == 1 ) 
                if ( abs(current_rad - abs(target_ew_1)) / current_rad <= 0.1 )
                    rad_ratio = 0.92; %0.85;
                else
                    rad_ratio = 0.95;
                end
            else
                angle_ew = mod(angle(all_ew_current_rad), 2*pi);
                idx      = find( abs(pi-angle_ew) > 2/180*pi, 1);
                if ( isempty(idx) )
                    rad_ratio = 0.95;
                else
                    rad_ratio = 1;
                end
            end
        
            rad = abs(target_ew_1) * ( 1 - 1.0e-3 );
    
            if ( abs(rad - current_rad) / rad < 1.0e-3 )
                rad = 0.95 * rad; 
            end
    
            if ( no_neg_real_ew <= no_pos_real_ew )
        
                if ( abs(min_rad_pos_real_ew - rad) >= (current_rad - rad)*0.45 )
                    max_wanted_ew = max(10, min(30, no_pos_real_ew));
            
                    if ( abs(min_rad_pos_real_ew) / rad >= 1.5 )
                        rad = ( rad + 3*abs(min_rad_pos_real_ew) ) / 4; 
                    end
        
                else
                    rad           = rad * rad_ratio;
                    max_wanted_ew = 10; 
                end
            else  
                     
                if ( abs(min_rad_neg_real_ew - rad) >= abs(current_rad - rad)*0.45 )
                    max_wanted_ew = max(15, min(32, no_neg_real_ew));
            
                    if ( abs(min_rad_neg_real_ew) / rad >= 1.5 )
                        rad       = ( 3 * rad + abs(min_rad_neg_real_ew) ) / 4; 
                    elseif ( abs(min_rad_neg_real_ew - rad)/rad < 1.0e-3 )
                        rad       = rad * rad_ratio; 
                    end
             
                else
                    rad           = rad * rad_ratio;
                    max_wanted_ew = 10; 
                end
            end
        
            current_diff_rad = current_rad - rad;
            fprintf('current_diff_rad = %11.4e,  diff_rad = %11.4e \n',current_diff_rad, diff_rad);
            if ( abs(current_diff_rad) > 0.2 )
                %
                % Modify at 2017/1/11
                %
                idx_radius_moving =  radius_moving > 0 ;
                rad               = max(rad, radius_moving(idx_radius_moving));
                %rad               = current_rad - 0.2;
            elseif ( current_diff_rad / diff_rad < 0.05 )
                rad = rad * 0.9;
                fprintf('Reduce the radius from %11.4e to %11.4e \n', rad/0.9, rad);
            end
            
            idx = find( abs(prediction_pt) > 0 );
            if ( ~isempty(idx) )
                prediction_radius = min(abs(prediction_pt(idx)));
                fprintf('prediction_radius = %11.4e, mean = %11.4e \n', prediction_radius, mean(abs(prediction_pt(idx))));
                
                if ( rad / prediction_radius < 1.25 || length(sort_idx) == 1 )
                    if ( rad / prediction_radius < 1.25 )
                        rad = min( rad, prediction_radius );
                    else
                        rad = 0.6 * rad + 0.4 * prediction_radius;
                    end
                end
            end
            
        case 3
            if ( current_rad >= 1.0e-3 )
                rad       = mean(abs(ew_sort2(1:end-1,1))) * 0.985;
            else
                rad       = mean(abs(ew_sort2(1:end-1,1))) * 1.05;
                if ( rad >= current_rad )
                    rad       = mean(abs(ew_sort2(1:end-1,1))) * 0.99;
                end
            end
            
            max_wanted_ew = 40;
            fprintf('classing case with radius = %11.4e \n', rad);
        
    end

end

idxx           = find(abs(ew_inside_c_new) <= current_rad+diff_rad/3);
    if ( ~isempty(idxx) ) 
        no_ew_in       = length(idxx);
        target_ew      = ew_inside_c_new(idxx);
        ang_ew         = mod(angle(target_ew),2*pi); 
        [~,idx_sort]   = sort(ang_ew);
        target_ew      = target_ew(idx_sort);
        ang_ew         = ang_ew(idx_sort);
        diff_ang       = abs(ang_ew(2:end) - ang_ew(1:end-1));
        idx_diff_ang   = find(diff_ang > 15/180*pi);
        rat_ew         = 0;
        
        if ( ~isempty(idx_diff_ang) )
            if ( length(idx_diff_ang) == 1 )
                no_target_ew            = idx_diff_ang;
            else
                no_target_ew                 = zeros(3,1);
                no_target_ew(1,1)            = idx_diff_ang(1);
                [no_target_ew(2,1), max_idx] = max(idx_diff_ang(2:end)-idx_diff_ang(1:end-1));
                no_target_ew(3,1)            = length(target_ew) - idx_diff_ang(end);
            end
            [no_cluster_ew, idx2]        = max(no_target_ew);
    
            switch idx2
                case 1
                    target_rad = min(abs(target_ew(1:idx_diff_ang(1))));
                case 2
                    if ( max_idx > 1 )
                        target_rad = min(abs(target_ew(idx_diff_ang(max_idx-1):idx_diff_ang(max_idx))));
                    else
                        target_rad = min(abs(target_ew(idx_diff_ang(max_idx):idx_diff_ang(max_idx+1))));
                    end
                case 3
                    target_rad = min(abs(target_ew(idx_diff_ang(end)+1:end)));
            end
     
            rat_ew         = no_cluster_ew/no_ew_in;
    
            fprintf('no_cluster_ew = %3.0f, ratio_ew = %11.4e in determine_radius_circle\n', no_cluster_ew, rat_ew);
            %target_rad
    
%             if ( rat_ew > 0.55 && no_cluster_ew >= 40) 
%                 if ( (current_rad - target_rad) / target_rad > 0.3 )
%                     rad = 1.02 * target_rad;
%                     fprintf('rat_ew = %11.4e, no_cluster_ew = %3.0f in determine_radius_circle (i)\n', rat_ew, no_cluster_ew);
%                 else
%                     rad = 1.005 * target_rad; %min(abs(target_ew(idx_cluster_ew))); %0.995 * min(abs(target_ew(idx_cluster_ew))); 
%                     fprintf('rat_ew = %11.4e, no_cluster_ew = %3.0f in determine_radius_circle (ii)\n', rat_ew, no_cluster_ew);
%                 end
%             elseif ( no_cluster_ew >= 40 ) 
%                 rad = target_rad;
%             end
        end
    end
    
%fprintf('current_rad = %11.4e in determine_radius_circle \n', current_rad);
if ( current_rad >= 1.5e-3 ) %( current_rad >= 8.0e-4 ) 
    if ( (current_rad - rad) / diff_rad < 1.4 )
        rad = rad * 0.9;
    end
    if ( rat_ew > 0.4 && (target_rad - rad)/rad > 0.2 && no_cluster_ew > 20)
        rad = target_rad;
    end
else
%     idxx           = find(abs(ew_inside_c_new) <= current_rad);
    if ( ~isempty(idxx) ) 
%         no_ew_in       = length(idxx);
%         target_ew      = ew_inside_c_new(idxx);
%         ang_ew         = mod(angle(target_ew),2*pi); 
%         [~,idx_sort]   = sort(ang_ew);
%         target_ew      = target_ew(idx_sort);
%         ang_ew         = ang_ew(idx_sort);
%         diff_ang       = abs(ang_ew(2:end) - ang_ew(1:end-1));
%         idx_diff_ang   = find(diff_ang > 15/180*pi);
%     
        if ( ~isempty(idx_diff_ang) )
%             if ( length(idx_diff_ang) == 1 )
%                 no_target_ew            = idx_diff_ang;
%             else
%                 no_target_ew                 = zeros(3,1);
%                 no_target_ew(1,1)            = idx_diff_ang(1);
%                 [no_target_ew(2,1), max_idx] = max(idx_diff_ang(2:end)-idx_diff_ang(1:end-1));
%                 no_target_ew(3,1)            = length(target_ew) - idx_diff_ang(end);
%             end
%             [no_cluster_ew, idx2]        = max(no_target_ew);
%     
%             switch idx2
%                 case 1
%                     target_rad = min(abs(target_ew(1:idx_diff_ang(1))));
%                 case 2
%                     target_rad = min(abs(target_ew(idx_diff_ang(max_idx-1):idx_diff_ang(max_idx))));
%                 case 3
%                     target_rad = min(abs(target_ew(idx_diff_ang(end)+1:end)));
%             end
%      
%             rat_ew         = no_cluster_ew/no_ew_in;
%     
%             %fprintf('ratio_ew = %11.4e in determine_radius_circle\n', rat_ew);
    
            if ( rat_ew > 0.55 && no_cluster_ew >= 40) 
                if ( (current_rad - target_rad) / target_rad > 0.3 )
                    rad = 1.02 * target_rad;
                    fprintf('rat_ew = %11.4e, no_cluster_ew = %3.0f in determine_radius_circle (i)\n', rat_ew, no_cluster_ew);
                else
                    rad = 1.005 * target_rad; %min(abs(target_ew(idx_cluster_ew))); %0.995 * min(abs(target_ew(idx_cluster_ew))); 
                    fprintf('rat_ew = %11.4e, no_cluster_ew = %3.0f in determine_radius_circle (ii)\n', rat_ew, no_cluster_ew);
                end
            elseif ( no_cluster_ew >= 40 ) 
                rad = target_rad;
            end
        end
    end
end

end

