
clear

Eng_eps = 3;

switch Eng_eps 
        case 0
            load ../0804/EW_data_0
            plot_eig = 'yes';
            
        case 0.5
            load ../0804/EW_data_0p5
            plot_eig = 'yes';
            
        case 1
            load ../0804/EW_data_1p0
            plot_eig = 'yes';
            
        case 1.5
            load ../0804/EW_data
            plot_eig = 'yes';
            
        case 2
            load ../0804/EW_data_2p0
            plot_eig = 'yes';
            
        case 2.5
            load ../0804/EW_data_2p5
            plot_eig = 'yes';
            
        case 3
            load ../0804/EW_data_3p0
            plot_eig = 'yes';
            
        case 3.5
            load ../0804/EW_data_3p5
            plot_eig = 'yes';
            
        case 4
            load ../0804/EW_data_4p0
            plot_eig = 'yes';
            
        case 4.5
            load ../0804/EW_data_4p5
            plot_eig = 'yes';
            
        case 5
            load ../0804/EW_data_5p0
            plot_eig = 'yes';
    
        case 5.5
            load ../0804/EW_data_55
            plot_eig = 'yes';
            
        case 6
            load ../0804/EW_data_60
            plot_eig = 'yes';
            
        case 7
            load ../0804/EW_data_7p0
            plot_eig = 'yes';
            
        case -2
            load ../v1_0919/EW_data_n2p0
            plot_eig = 'yes';
            
        case -5.1
            load ../v2_0924/EW_data_n5p1
            plot_eig = 'yes';
            
        case -6.5
            load ../v1_0919/EW_data_n6p5
            plot_eig = 'yes';
            
        otherwise
            plot_eig = 'no';
            
end

switch Eng_eps 
        case 0
            load tmp_dat_eps0 
            
        case 0.5
            load tmp_dat_eps0p5  
            
        case 1
            load tmp_dat_eps1p0  
            
        case 1.5
            load tmp_dat_eps1p5  %tmp_dat_eps1p5_23
            
        case 2
            load tmp_dat_eps2p0
            
        case 2.5
            load tmp_dat_eps2p5 
            
        case 3
            load tmp_dat_eps3p0 
            
        case 3.5
            load tmp_dat_eps3p5 
            
        case 4
            load tmp_dat_eps4p0 
            
        case 4.5
            load tmp_dat_eps4p5 
            
        case 5
            load tmp_dat_eps5p0 
    
        case 5.5
            load tmp_dat_eps5p5 
            
        case 6
            load tmp_dat_eps6p0 
            
        case 7
            load tmp_dat_eps7p0 
            
        case -2
            load tmp_dat_epsn2p0
            
        case -5.1
            load tmp_dat_epsn5p1
            
        case -6.5
            load tmp_dat_epsn6p5
            
end
    
idx_target_ew = find(abs(ew)<= 1+1.0e-9 & imag(ew) >= -1.0e-9);
target_ew     = ew(idx_target_ew);
[~,idx_ew]    = sort(real(target_ew));

no_ew_eig             = length(target_ew);
    idx_unit_ew_eig       = find(abs(abs(target_ew)-1) < 1.0e-10);
    idx_cmp_ew_eig        = find(abs(imag(target_ew)) > 1.0e-9 & abs(abs(target_ew)-1) >= 1.0e-10);
    idx_real_ew_eig       = setdiff(1:no_ew_eig, [idx_cmp_ew_eig; idx_unit_ew_eig]);
    no_target_ew_on_c_eig = length(idx_unit_ew_eig);
    no_target_cmp_ew_eig  = length(idx_cmp_ew_eig);
    no_target_real_ew_eig = length(idx_real_ew_eig);
total_ew_eig         = no_target_ew_on_c_eig + 2*no_target_cmp_ew_eig + no_target_real_ew_eig;

no_ew             = length(ew_inside_c);
tt                = [ew_unit_c(:,1); ew_unit_c(:,2)];
idx_unit_ew       = find( imag(tt) >= -1.0e-7 );
    idx_unit_ew2      = find(abs(abs(ew_inside_c)-1) < 1.0e-10);
    idx_cmp_ew        = find(abs(imag(ew_inside_c)) > 1.0e-9 & abs(abs(ew_inside_c)-1) >= 1.0e-10);
    idx_real_ew       = setdiff(1:no_ew, [idx_cmp_ew; idx_unit_ew2]);
    no_target_ew_on_c = length(idx_unit_ew);
    no_target_cmp_ew  = length(idx_cmp_ew);
    no_target_real_ew = length(idx_real_ew);
total_ew         = no_target_ew_on_c + 2*no_target_cmp_ew + no_target_real_ew;

error_un_ew   = no_target_ew_on_c_eig - no_target_ew_on_c;
error_cmp_ew  = no_target_cmp_ew_eig - no_target_cmp_ew;
error_real_ew = no_target_real_ew_eig - no_target_real_ew;

no_loss_ew = error_un_ew + error_cmp_ew + error_real_ew;
if ( no_loss_ew == 0 )
    fprintf('All eigenvalues are computed \n');
else
    loss_ew = zeros(no_loss_ew,1);
    if ( error_un_ew > 0 )
        [~,idx1] = sort(real(target_ew(idx_unit_ew_eig)));
        [~,idx2] = sort(real(tt(idx_unit_ew)));
        tmp_ew   = target_ew(idx_unit_ew_eig(idx1));
        for ii = 1:error_un_ew
            error1   = tmp_ew(1:no_target_ew_on_c) - tt(idx_unit_ew(idx2));
            idx3     = find(abs(error1) > 1.0e-7);
            
            loss_ew(ii)           = tmp_ew(idx3(1));
            tmp_ew(idx3(1):end-1) = tmp_ew(idx3(1)+1:end);
        end 
    end
    
    if ( error_real_ew > 0 )
        [~,idx1] = sort(real(target_ew(idx_real_ew_eig)));
        [~,idx2] = sort(real(ew_inside_c(idx_real_ew)));
        tmp_ew   = target_ew(idx_real_ew_eig(idx1));
        %error1   = abs(target_ew(idx_real_ew_eig(idx1(1:no_target_real_ew))) - ew_inside_c(idx_real_ew(idx2)));
        
        for ii = 1:error_real_ew
            error1   = abs(tmp_ew(1:no_target_real_ew) - ew_inside_c(idx_real_ew(idx2)));
            error1   = error1./abs(tmp_ew(1:no_target_real_ew));
            idx3     = find(error1 > 1.0e-7);
        
            loss_ew(error_un_ew+ii) = tmp_ew(idx3(1));
            tmp_ew(idx3(1):end-1)   = tmp_ew(idx3(1)+1:end);
        end
    end
    
    if ( error_cmp_ew > 0 ) 
        tmp_eig     = target_ew(idx_cmp_ew_eig);
        tmp_ew      = ew_inside_c(idx_cmp_ew); 
        [~,idx1]    = sort(real(tmp_eig));
        [~,idx2]    = sort(real(tmp_ew));
        tmp_eig     = tmp_eig(idx1);
        tmp_ew      = tmp_ew(idx2);
        test        = tmp_eig;
        
        idx_loss_ew = zeros(no_target_cmp_ew_eig,1);
        for ii = 1:no_target_cmp_ew
            tt                  = abs(tmp_ew(ii) - tmp_eig) / abs(tmp_ew(ii));
            [~,idx2]            = min(tt);
            idx_loss_ew(idx2,1) = 1;
            [~,idxtest]         = min(abs(tmp_ew(ii) - test) / abs(tmp_ew(ii)));
            if ( idx2 ~= idxtest )
                ii
                tmp_ew(ii)
                tmp_eig(idx2)
                test(idxtest)
            end
            test(idx2,1)        = -10;
        end
        idx_loss_target_ew = find( idx_loss_ew == 0);
        
        if ( length(idx_loss_target_ew) ~= error_cmp_ew )
            fprintf('Error for finding complex eigenvalues \n');
        end
        loss_ew(error_un_ew+error_real_ew+1:error_un_ew+error_real_ew+length(idx_loss_target_ew)) = tmp_eig(idx_loss_target_ew);
        
%         for ii = 1:error_cmp_ew
%             error1   = abs(tmp_ew(1:no_target_cmp_ew) - ew_inside_c(idx_cmp_ew(idx2)));
%             error1   = error1./abs(tmp_ew(1:no_target_cmp_ew));
%             idx3     = find(error1 > 1.0e-7);
%         
%             loss_ew(error_un_ew+error_real_ew+ii) = tmp_ew(idx3(1));
%             tmp_ew(idx3(1):end-1)                 = tmp_ew(idx3(1)+1:end);
%         end
    end
end