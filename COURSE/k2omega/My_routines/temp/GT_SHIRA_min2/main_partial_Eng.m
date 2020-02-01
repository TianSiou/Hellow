clear
format long e
%clc
% delete rslt.txt
% diary rslt.txt
% diary on
global ew

Eng_eps    = -5.1; %-2.0; %-2.0; %-6.5; %3.1; %0.5; %1.0; %0; %1.5; %7; %5.0; %4.5;

load_data  = 'no'; %'yes';
debug_flag = 'no';
max_circle = 28; %17; %15
plot_fig   = 'no';

load ../mat_datas/h00L
mtx_H0_large = A;
load ../mat_datas/h01L
mtx_H1_large = A;

load ../mat_datas/s00L
mtx_S0_large = A;
load ../mat_datas/s01L
mtx_S1_large = A;

fid           = fopen('parameter.txt','r');
Eng_eps_begin = fscanf(fid,'%g',1);
Eng_eps_end   = fscanf(fid,'%g',1);
Eng_eps_incre = fscanf(fid,'%g',1);
file_name     = fscanf(fid,'%s',1);
status        = fclose(fid);

fid           = fopen( file_name, 'w+');

for Eng_eps = Eng_eps_begin:Eng_eps_incre:Eng_eps_end
    fprintf('\n\n =============== Eng_eps = %11.4e ================== \n', Eng_eps);
dim_large = size(mtx_H0_large,1);
dim_n     = dim_large / 2;
mtx_H0    = mtx_H0_large(1:dim_n, 1:dim_n);
mtx_H1    = mtx_H1_large(1:dim_n,dim_n+1:2*dim_n);
mtx_S0    = mtx_S0_large(1:dim_n, 1:dim_n);
mtx_S1    = mtx_S1_large(1:dim_n,dim_n+1:2*dim_n);

mtx_A   = Eng_eps * mtx_S1 - mtx_H1;
mtx_Q   = Eng_eps * mtx_S0 - mtx_H0;

mtx_NME.mtx_A = Eng_eps * mtx_S1 - mtx_H1;
mtx_NME.mtx_Q = Eng_eps * mtx_S0 - mtx_H0;


% theta1           = pi/2 + 1.1e-3;
% theta2           = 1.5 * pi + 1.0e-3;
%theta            = linspace(0, pi, 8);
gamma            = 1; %exp(pi/12*1i); %[1; exp(1i*(theta(2:end-1)'+ 1.0e-3)); -1]; %1;%exp(1i*(theta(2:end-1)'+ 1.0e-3));
leng_r           = length(gamma);
CPUTime_LU       = zeros(leng_r,1);
CPUtime_JLanczos = zeros(leng_r,1);
CPUtime_gen_M    = zeros(leng_r,1);
CPUtime_eig      = zeros(leng_r,1);
NO_ew            = zeros(leng_r,1);
Es_ew_symplectic        = zeros(4*dim_n,leng_r);
ew_Hamiltonian   = zeros(2*dim_n,leng_r);

rad_circle       = zeros(max_circle,1);

if ( strcmp(plot_fig, 'yes') )
    switch Eng_eps 
        case 0 
            vidObj = VideoWriter('ew_circleBYcircle_0p0.mp4', 'MPEG-4');
            
        case 0.5
            vidObj = VideoWriter('ew_circleBYcircle_0p5.mp4', 'MPEG-4');
            
        case 1
            vidObj = VideoWriter('ew_circleBYcircle_1p0.mp4', 'MPEG-4');
            
        case 1.5
            vidObj = VideoWriter('ew_circleBYcircle_1p5.mp4', 'MPEG-4');
            
        case 2
            vidObj = VideoWriter('ew_circleBYcircle_2p0.mp4', 'MPEG-4');
            
        case 2.5
            vidObj = VideoWriter('ew_circleBYcircle_2p5.mp4', 'MPEG-4');
            
        case 3
            vidObj = VideoWriter('ew_circleBYcircle_3p0.mp4', 'MPEG-4');
            
        case 3.5
            vidObj = VideoWriter('ew_circleBYcircle_3p5.mp4', 'MPEG-4');
            
        case 4
            vidObj = VideoWriter('ew_circleBYcircle_4p0.mp4', 'MPEG-4');
            
        case 4.5
            vidObj = VideoWriter('ew_circleBYcircle_4p5.mp4', 'MPEG-4');
            
        case 5
            vidObj = VideoWriter('ew_circleBYcircle_5p0.mp4', 'MPEG-4');
    
        case 5.5
            vidObj = VideoWriter('ew_circleBYcircle_5p5.mp4', 'MPEG-4');
            
        case 6
            vidObj = VideoWriter('ew_circleBYcircle_6p0.mp4', 'MPEG-4');
            
        case 7
            vidObj = VideoWriter('ew_circleBYcircle_7p0.mp4', 'MPEG-4');
        
        case -6.5
            vidObj = VideoWriter('ew_circleBYcircle_n6p5.mp4', 'MPEG-4');
            %vidObj = VideoWriter('ew_circleBYcircle_n6p5.avi');
        
        case -2.0
            vidObj = VideoWriter('ew_circleBYcircle_n2p0.mp4', 'MPEG-4');
        
        case -5.1
            vidObj = VideoWriter('ew_circleBYcircle_n5p1_tmp.mp4', 'MPEG-4');
            
        otherwise
            vidObj = VideoWriter('ew_circleBYcircle.mp4', 'MPEG-4'); 
            
    end
    
    vidObj.FrameRate = 1;  % Default 30

    open(vidObj);
    
    theta = 0:0.001:2*pi;
    x     = cos(theta)+1i*sin(theta);

    figure(1) 
    plot(x,'r-')
    axis([-1.1, 1.1, -1.1, 1.1]);
    hold on

    figure(2) 
    plot(x,'r-')
    axis([-1.1, 1.1, -1.1, 1.1]);
    title(['\epsilon = ',num2str(Eng_eps)])
    hold on
    
else
    vidObj = [];
end

flag_full = 0;
if ( flag_full == 1 )
    mtx_M = [ full(mtx_A) zeros(dim_n); full(mtx_Q) -eye(dim_n)];
    mtx_L = [ zeros(dim_n) eye(dim_n); full(mtx_A') zeros(dim_n)];
    
    tic;
    [EV,EW] = eig(mtx_M, mtx_L, 'qz');
    CPUTime_qz = toc;
    fprintf('CPU time for QZ is %11.4e. \n\n', CPUTime_qz)
    
    ew = diag(EW);
    %save EW_data_0p5 ew

    Err = mtx_M * EV - mtx_L * EV * EW;
%     Res = norm(mtx_M * EV - mtx_L * EV * EW,1)
    res = zeros(2*dim_n,1);
    for ii = 1:2*dim_n
        res(ii,1) = norm(Err(:,ii));
    end
    
    save EW_data_n6p5 ew res
    
    plot_eig = 'yes';
    
else
    
    switch Eng_eps 
        case 0
            load ../mat_datas/EW_data_0
            plot_eig = 'yes';
            
        case 0.5
            load ../mat_datas/EW_data_0p5
            plot_eig = 'yes';
            
        case 1
            load ../mat_datas/EW_data_1p0
            plot_eig = 'yes';
            
        case 1.5
            load ../mat_datas/EW_data
            plot_eig = 'yes';
            
        case 2
            load ../mat_datas/EW_data_2p0
            plot_eig = 'yes';
            
        case 2.5
            load ../mat_datas/EW_data_2p5
            plot_eig = 'yes';
            
        case 3
            load ../mat_datas/EW_data_3p0
            plot_eig = 'yes';
            
        case 3.5
            load ../mat_datas/EW_data_3p5
            plot_eig = 'yes';
            
        case 4
            load ../mat_datas/EW_data_4p0
            plot_eig = 'yes';
            
        case 4.5
            load ../mat_datas/EW_data_4p5
            plot_eig = 'yes';
            
        case 5
            load ../mat_datas/EW_data_5p0
            plot_eig = 'yes';
    
        case 5.5
            load ../mat_datas/EW_data_55
            plot_eig = 'yes';
            
        case 6
            load ../mat_datas/EW_data_60
            plot_eig = 'yes';
            
        case 7
            load ../mat_datas/EW_data_7p0
            plot_eig = 'yes';
            
        case -2
            load ../mat_datas/EW_data_n2p0
            plot_eig = 'yes';
            
        case -5.1
            load ../mat_datas/EW_data_n5p1
            plot_eig = 'yes';
            
        case -6.5
            load ../mat_datas/EW_data_n6p5
            plot_eig = 'yes';
            
        otherwise
            plot_eig = 'no';
            
    end
end

if ( strcmp(plot_fig, 'yes') && strcmp(plot_eig, 'yes') )
    figure(1)
    plot(ew,'bo','LineWidth',1); 
    hold on

    figure(2)
    plot(ew,'bo','LineWidth',1); 
    hold on
end

[i_idx_A, j_idx_A, val_A] = find(mtx_NME.mtx_A);
[i_idx_Q, j_idx_Q, val_Q] = find(mtx_NME.mtx_Q);
            
nmz_A   = length(i_idx_A);
nmz_Q   = length(i_idx_Q);
i_idx_M = zeros(nmz_A+nmz_Q+dim_n,1);
j_idx_M = zeros(nmz_A+nmz_Q+dim_n,1);
val_M   = zeros(nmz_A+nmz_Q+dim_n,1);
            
i_idx_M(1:nmz_A,1) = i_idx_A;
j_idx_M(1:nmz_A,1) = j_idx_A;
val_M(1:nmz_A,1)   = val_A;
i_idx_M(1+nmz_A:nmz_A+nmz_Q,1) = dim_n + i_idx_Q;
j_idx_M(1+nmz_A:nmz_A+nmz_Q,1) = j_idx_Q;
val_M(1+nmz_A:nmz_A+nmz_Q,1)   = val_Q;
i_idx_M(1+nmz_A+nmz_Q:nmz_A+nmz_Q+dim_n,1) = (dim_n+1:2*dim_n)';
j_idx_M(1+nmz_A+nmz_Q:nmz_A+nmz_Q+dim_n,1) = (dim_n+1:2*dim_n)';
val_M(1+nmz_A+nmz_Q:nmz_A+nmz_Q+dim_n,1)   = -ones(dim_n,1);
            
mtx_M = sparse(i_idx_M, j_idx_M, val_M);

i_idx_L = zeros(nmz_A+dim_n,1);
j_idx_L = zeros(nmz_A+dim_n,1);
val_L   = zeros(nmz_A+dim_n,1);
            
i_idx_L(1:dim_n,1) = (1:dim_n)';
j_idx_L(1:dim_n,1) = (dim_n+1:2*dim_n)';
val_L(1:dim_n,1)   = ones(dim_n,1);
i_idx_L(dim_n+1:dim_n+nmz_A,1) = dim_n+j_idx_A;
j_idx_L(dim_n+1:dim_n+nmz_A,1) = i_idx_A;
val_L(dim_n+1:dim_n+nmz_A,1)   = val_A;
     
mtx_L     = sparse(i_idx_L, j_idx_L, val_L);
            
deflation = 'yes';
    
if ( strcmp(load_data,'no') )
%     switch deflation 
%         case 'yes'
            
            rad                = 1;
            rad_old            = 1;
            max_wanted_ew      = 45;
            radius_case        = 'normal';
            no_circle          = 1;
            weighting_choosing = 'right'; 

            ew_unit_c        = []; 
            ev_on_c          = []; 
            ew_inside_c      = [];
            ev_inside_c      = []; 
            ev_outside_c     = []; 
            ref_ew_inside_c  = []; 
            ref_ev_inside_c  = [];
            
            [ ew_unit_c, ev_on_c, ew_inside_c, ev_inside_c, j_no_deuplicate_ew, ...
                rad, ev_outside_c, ref_ew_inside_c, ref_ev_inside_c, max_wanted_ew, flag_conv, ...
                ew_inner, radius_case, weighting_choosing, jj ] = compute_1st_circle_ew( dim_n, mtx_NME, ...
                mtx_M, mtx_L, rad, rad_old, max_wanted_ew, radius_case, vidObj, weighting_choosing, ...
                ew_unit_c, ev_on_c, ew_inside_c, ev_outside_c, ref_ew_inside_c, ref_ev_inside_c, plot_fig, deflation );

%             tmp_vec2 = zeros(size(ew_unit_c,1),2);
%             for j2 = 1:size(ew_unit_c,1)
%                 for i2 = 1:2
%                     tmp_vec2(j2,i2) = 1i * (ev_on_c(1:dim_n,j2,i2)'*(2 * ew_unit_c(j2,i2) * ...
%                         (mtx_NME.mtx_A.'*ev_on_c(1:dim_n,j2,i2)) - mtx_NME.mtx_Q*ev_on_c(1:dim_n,j2,i2)));
% %                     fprintf('d(%24.16e+(%24.16e)*1i,%1.0f) = %11.4e+(%11.4e)*1i \n',real(ew_unit_c(j2,i2)), ...
% %                         imag(ew_unit_c(j2,i2)),i2, real(tmp_vec2(j2,i2)), imag(tmp_vec2(j2,i2)))
%                 end
%             end
%         
%             idx = find( abs(imag(tmp_vec2)) > 1.0e-10 );
%             if ( ~isempty(idx) )
%                 fprintf('Error in choosing unimodular eigenvalues with index \n');
%                 display(idx)
%             else
%                 idx1           = find(real(tmp_vec2(:,1)) > 0 & real(tmp_vec2(:,2)) < 0);
%                 idx2           = find(real(tmp_vec2(:,2)) > 0 & real(tmp_vec2(:,1)) < 0);
%                 if ( length(idx1) + length(idx2) ~= size(ew_unit_c,1))
%                     fprintf('Error in choosing unimodular eigenvalues \n');
%                 else
%                     target_ev_on_c = zeros(2*dim_n,size(ew_unit_c,1));
%                     target_unit_ew = [ ew_unit_c(idx1,1); ew_unit_c(idx2,2) ];
%                     target_ev_on_c(:,1             :length(idx1)             ) = ev_on_c(:,idx1,1);
%                     target_ev_on_c(:,1+length(idx1):length(idx1)+length(idx2)) = ev_on_c(:,idx2,2);
%                     %target_ev_on_c = [ ev_on_c(:,idx1,1); ev_on_c(:,idx2,2) ];
%                 end
%             end
        
%         case 'no'
%             [ ew_unit_c, ev_on_c, ew_inside_c, ev_inside_c, no_deuplicate_ew ] = compute_ew_near_unit( dim_n, mtx_NME );
%     end

    switch Eng_eps 
        case 0
            save near_c_data_eps_0  rad ew_unit_c ev_on_c ew_inside_c ev_inside_c ev_outside_c ew_inner
            
        case 0.5
            save near_c_data_eps_0p5   rad ew_unit_c ev_on_c ew_inside_c ev_inside_c ev_outside_c ew_inner
            
        case 1
            save near_c_data_eps_1p0   rad ew_unit_c ev_on_c ew_inside_c ev_inside_c ev_outside_c ew_inner
            
        case 1.5
            save near_c_data_eps_1p5  rad ew_unit_c ev_on_c ew_inside_c ev_inside_c ev_outside_c ew_inner
            
        case 2
            save near_c_data_eps_2p0  rad ew_unit_c ev_on_c ew_inside_c ev_inside_c ev_outside_c ew_inner
            
        case 2.5
            save near_c_data_eps_2p5  rad ew_unit_c ev_on_c ew_inside_c ev_inside_c ev_outside_c ew_inner
            
        case 3
            save near_c_data_eps_3p0  rad ew_unit_c ev_on_c ew_inside_c ev_inside_c ev_outside_c ew_inner
            
        case 3.5
            save near_c_data_eps_3p5  rad exceptional rad_check ew_unit_c ev_on_c ew_inside_c ev_inside_c ev_outside_c ew_inner
            
        case 4
            save near_c_data_eps_4p0  rad ew_unit_c ev_on_c ew_inside_c ev_inside_c ev_outside_c ew_inner
            
        case 4.5
            save near_c_data_eps_4p5  rad  ew_unit_c ev_on_c ew_inside_c ev_inside_c ev_outside_c ew_inner
            
        case 5
            save near_c_data_eps_5p0  rad ew_unit_c ev_on_c ew_inside_c ev_inside_c ev_outside_c ew_inner
    
        case 5.5
            save near_c_data_eps_5p5  rad ew_unit_c ev_on_c ew_inside_c ev_inside_c ev_outside_c ew_inner
            
        case 6
            save near_c_data_eps_6p0  rad ew_unit_c ev_on_c ew_inside_c ev_inside_c ev_outside_c ew_inner
            
        case 7
            save near_c_data_eps_7p0  rad ew_unit_c ev_on_c ew_inside_c ev_inside_c ev_outside_c ew_inner
            
        case -2
            save near_c_data_eps_n2p0  rad ew_unit_c ev_on_c ew_inside_c ev_inside_c ev_outside_c ew_inner
            
        case -5.1
            save near_c_data_eps_n5p1  rad ew_unit_c ev_on_c ew_inside_c ev_inside_c ev_outside_c ew_inner
            
        case -6.5
            save near_c_data_eps_n6p5  rad ew_unit_c ev_on_c ew_inside_c ev_inside_c ev_outside_c ew_inner
            
    end
    
    %save near_c_data rad exceptional rad_check ew_unit_c ev_on_c ew_inside_c ev_inside_c ev_outside_c

else
    
    switch Eng_eps 
        case 0
            load near_c_data_eps_0 
            
        case 0.5
            load near_c_data_eps_0p5  
            
        case 1
            load near_c_data_eps_1p0  
            
        case 1.5
            load near_c_data_eps_1p5 
            
        case 2
            load near_c_data_eps_2p0 
            
        case 2.5
            load near_c_data_eps_2p5 
            
        case 3
            load near_c_data_eps_3p0 
            
        case 3.5
            load near_c_data_eps_3p5 
            
        case 4
            load near_c_data_eps_4p0 
            
        case 4.5
            load near_c_data_eps_4p5 
            
        case 5
            load near_c_data_eps_5p0 
    
        case 5.5
            load near_c_data_eps_5p5 
            
        case 6
            load near_c_data_eps_6p0 
            
        case 7
            load near_c_data_eps_7p0 
            
        case -2
            load near_c_data_eps_n2p0 
            
        case -5.1
            load near_c_data_eps_n5p1
            
        case -6.5
            load near_c_data_eps_n6p5 
            
    end
    
%     load near_c_data_eps_0
    if (strcmp(plot_fig, 'yes') )    
        figure(2)
        plot(ew_unit_c(:,1),'g+','LineWidth',1); 
        plot(ew_unit_c(:,2),'g+','LineWidth',1);
        plot(ew_inside_c(:,1),'g+','LineWidth',1);
    end
end

if ( size(ew_unit_c,1) > 0 )
    %
    % Find the target eigenvalues on the unit circle
    %
            tmp_vec2 = zeros(size(ew_unit_c,1),2);
            for j2 = 1:size(ew_unit_c,1)
                for i2 = 1:2
                    tmp_vec2(j2,i2) = 1i * (ev_on_c(1:dim_n,j2,i2)'*(2 * ew_unit_c(j2,i2) * ...
                        (mtx_NME.mtx_A.'*ev_on_c(1:dim_n,j2,i2)) - mtx_NME.mtx_Q*ev_on_c(1:dim_n,j2,i2)));
%                     fprintf('d(%24.16e+(%24.16e)*1i,%1.0f) = %11.4e+(%11.4e)*1i \n',real(ew_unit_c(j2,i2)), ...
%                         imag(ew_unit_c(j2,i2)),i2, real(tmp_vec2(j2,i2)), imag(tmp_vec2(j2,i2)))
                end
            end
        
            idx = find( abs(imag(tmp_vec2)) > 1.0e-5 );
            if ( ~isempty(idx) )
                fprintf('Error in choosing unimodular eigenvalues with index \n');
                display(idx)
            else
                idx1           = find(real(tmp_vec2(:,1)) > 0 & real(tmp_vec2(:,2)) < 0);
                idx2           = find(real(tmp_vec2(:,2)) > 0 & real(tmp_vec2(:,1)) < 0);
                if ( length(idx1) + length(idx2) ~= size(ew_unit_c,1))
                    fprintf('Error in choosing unimodular eigenvalues \n');
                else
                    target_ev_on_c = zeros(2*dim_n,size(ew_unit_c,1));
                    target_unit_ew = [ ew_unit_c(idx1,1); ew_unit_c(idx2,2) ];
                    target_ev_on_c(:,1             :length(idx1)             ) = ev_on_c(:,idx1,1);
                    target_ev_on_c(:,1+length(idx1):length(idx1)+length(idx2)) = ev_on_c(:,idx2,2);
                    %target_ev_on_c = [ ev_on_c(:,idx1,1); ev_on_c(:,idx2,2) ];
                end
            end
    %
    % =================================
    %
    
    %n           = size(ev_inside_c,1);
    no_ew             = length(ew_inside_c);
    idx_unit_ew       = find(abs(abs(ew_inside_c)-1) < 1.0e-10);
    idx_cmp_ew        = find(abs(imag(ew_inside_c)) > 1.0e-9 & abs(abs(ew_inside_c)-1) >= 1.0e-10);
    idx_real_ew       = setdiff(1:no_ew, [idx_cmp_ew; idx_unit_ew]);
    no_target_ew_on_c = size(ew_unit_c,1); %length(idx_unit_ew);
    no_target_cmp_ew  = length(idx_cmp_ew);
    no_target_real_ew = length(idx_real_ew);

    %no_target_ew_on_c = length(idx1) + length(idx2);
    idx_in_upp_c      = find(imag(ew_inside_c) >= 0);
    no_target_ew_in_c = length(idx_in_upp_c);

    fprintf('---- begin to compute inner eigenvalues no_unit_ew = %3.0f and no_ew_in_c = %4.0f ------ \n', ...
        no_target_ew_on_c, no_target_cmp_ew*2+no_target_real_ew);

    ev_all_outside_c                    = zeros(dim_large, dim_n);
    no_outside_ev                       = size(ev_outside_c,2);
    ref_ev_outside_c                    = ev_outside_c;
    ref_ew_inside_c                     = ew_inside_c;
    ref_ev_inside_c                     = ev_inside_c;
    ev_all_outside_c(:,1:no_outside_ev) = ev_outside_c;

    idx_neg_real_ew  = find(abs(imag(ew_inside_c)) <= 1.0e-8 & angle(ew_inside_c) >= pi - 3.0e-2 ...
        & abs(ew_inside_c) >= rad);

    rad_circle(1,1)  = 1;
    rad_old          = 1; 

    if ( isempty(idx_neg_real_ew) )
        min_rad_neg_real_ew = rad_old;
    else
        min_rad_neg_real_ew = max(real(ew_inside_c(idx_neg_real_ew)));
    end

    if ( length(ew_inner) > 1 )
        rad             = max(rad, 0.65);
        if ( max_circle > 1 )
            rad_circle(2,1) = max(rad_circle(2,1), 0.65);
        end
        max_wanted_ew   = 15;
        m2              = 8;
        radius_case     = 'normal';
    else
        rad_circle(2,1) = 0.7 * rad_circle(2,1) + 0.3 * rad_circle(1,1);
        rad             = 0.7 * rad + 0.3 * rad_old;
        max_wanted_ew   = 45;
        m2              = max_wanted_ew;
        radius_case     = 'weighting';
    end


if ( strcmp(debug_flag,'yes') )
    
    weighting_choosing = 'left';
    
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
% load tmp_dat_eps0p5_20

    if ( strcmp(plot_fig, 'yes') )
        figure(2)
        plot(ew_inside_c(:,1),'g+','LineWidth',1);
        axis([ -rad_old rad_old -rad_old rad_old])

        figure(1)
        axis([ -rad_old rad_old -rad_old rad_old])
    end
    
    no_circle          = no_circle + 1;
    %weighting_choosing = 'left';
else
    no_circle          = 1;
    weighting_choosing = 'right';
end

no_circle_old = no_circle;
flag_stop     = 1;
no_cong_ew    = no_target_ew_on_c + 2 * no_target_cmp_ew + no_target_real_ew;

while ( no_circle <= max_circle && flag_stop == 1 )
%for no_circle = be_circle:30 
    
    if ( strcmp(plot_fig, 'yes') )
        figure(2) 
        plot(rad*x,'r-'); 
        tt = num2str(no_circle);
        text(rad, 0, tt)
        if ( strcmp(debug_flag,'yes') && no_circle <= no_circle_old+1 )
            axis([ -rad_old rad_old -rad_old rad_old])
        elseif ( no_circle == 1 )
            axis([ -rad_old rad_old -rad_old rad_old])
        elseif ( no_circle == 2 && abs(Eng_eps - 5.5) <= 1.0e-10 ) 
            axis([ -rad_circle(no_circle+1,1)-0.03 rad_circle(no_circle+1,1)+0.03 -rad_circle(no_circle+1,1)-0.03 rad_circle(no_circle+1,1)+0.03])
        elseif ( no_circle == 3 && abs(Eng_eps - 5.5) <= 1.0e-10 ) 
            axis([ -rad_circle(no_circle-1,1)+0.1 rad_circle(no_circle-1,1)-0.1 -rad_circle(no_circle-1,1)+0.1 rad_circle(no_circle-1,1)-0.1])
        elseif ( no_circle == 3 )  
            axis([ -rad_circle(no_circle-2,1) rad_circle(no_circle-2,1) -rad_circle(no_circle-2,1) rad_circle(no_circle-2,1)])
        else
            axis([ -rad_circle(no_circle-1,1) rad_circle(no_circle-1,1) -rad_circle(no_circle-1,1) rad_circle(no_circle-1,1)])
        end
        drawnow
        hold on
    end

    fprintf('\n ======= Begin %2.0f-th interior circle with no_unit_ew = %3.0f and no_ew_in_c = %4.0f ========= \n', no_circle, ...
        no_target_ew_on_c, no_target_cmp_ew*2+no_target_real_ew);
    
     [ ew_unit_c, ev_on_c, ew_inside_c_new, ev_inside_c_new, j_no_deuplicate_ew,  ...
        rad_new, ref_ev_outside_c, ref_ew_inside_c, ref_ev_inside_c, max_wanted_ew, flag_conv, ew_inner, ...
        radius_case, weighting_choosing, jj, no_target_ew_on_c, no_target_cmp_ew, no_target_real_ew ] = ...
        compute_jth_circle_ew_1228( dim_n, mtx_NME, mtx_M, mtx_L, rad, rad_old, max_wanted_ew, radius_case, ...
        vidObj, weighting_choosing, ew_unit_c, ev_on_c, ew_inside_c, ev_inside_c, ref_ev_outside_c, ...
        ref_ew_inside_c, ref_ev_inside_c, plot_fig, no_target_ew_on_c, no_target_cmp_ew, no_target_real_ew, deflation);
     
     prev_no_cong_ew = no_cong_ew;
     no_cong_ew      = no_target_ew_on_c + 2 * no_target_cmp_ew + no_target_real_ew;
     if ( no_cong_ew == dim_n || (no_cong_ew == prev_no_cong_ew && rad < 1.0e-3))
         flag_stop = 0;
     end
     
     if ( flag_conv == 0 )
         ew_inside_c               = ew_inside_c_new;
         ev_inside_c               = ev_inside_c_new;
         rad_old                   = rad;
         rad                       = rad_new;
         rad_circle(no_circle+2,1) = rad_new;
     
         idx_in_upp_c              = find(imag(ew_inside_c) >= 0);
         no_target_ew_in_c         = length(idx_in_upp_c);
         
         no_out_ew                                                   = size(ref_ev_outside_c,2);
         ev_all_outside_c(:,1+no_outside_ev:no_outside_ev+no_out_ew) = ref_ev_outside_c;
         no_outside_ev                                               = no_outside_ev+no_out_ew;
     
     end
     
     no_circle = no_circle + 1;
     
end

flag1             = 1;

[ X_sol, RelErr ] = form_eigspace( mtx_NME, mtx_M, mtx_L, ew_inside_c, ev_inside_c, ev_all_outside_c, ew_unit_c, ev_on_c );

no_loss_ew = dim_n - (no_target_ew_on_c + no_target_cmp_ew*2 + no_target_real_ew);
fprintf('Eng_eps(kk,1) = %11.4e; RelErr(kk,1) = %11.4e; no_loss_ew(kk) = %3.0f; kk = kk + 1; \n', Eng_eps, RelErr, no_loss_ew);
fprintf(fid,'Eng_eps(kk,1) = %11.4e; RelErr(kk,1) = %11.4e; no_loss_ew(kk) = %3.0f; kk = kk + 1; \n', Eng_eps, RelErr, no_loss_ew);

% if ( no_cong_ew ~= dim_n && flag1 == 0) 
%     if ( no_cong_ew < dim_n )
%         %U_eig_space = form_eigspace( ew_inside_c, ev_inside_c, ev_all_outside_c );
%         
%         no_ew       = length(ew_inside_c);
%         idx_unit_ew = find(abs(abs(ew_inside_c)-1) < 1.0e-10);
%         idx_cmp_ew  = find(abs(imag(ew_inside_c)) > 1.0e-9 & abs(abs(ew_inside_c)-1) >= 1.0e-10);
%         idx_real_ew = setdiff(1:no_ew, [idx_cmp_ew; idx_unit_ew]);
%         mm_real     = length(idx_real_ew);
%         mm_cmp      = length(idx_cmp_ew);
%         mm_unit_ew  = length(idx_unit_ew);
%         mm          = mm_real + 2 * mm_cmp;
%         U_in        = zeros(dim_large,mm);
%         U_out       = zeros(dim_large,mm);
%         
%         U_in (:,        1:mm_real) = real(ev_inside_c(:,idx_real_ew));
%         U_out(:,        1:mm_real) = real(ev_all_outside_c(:,idx_real_ew));
%         for ii = 1:mm_real
%             U_in(:,ii) = U_in(:,ii) / norm(U_in(:,ii));
%             U_out(:,ii) = U_out(:,ii) / (U_in(:,ii)' * [U_out(dim_n+1:dim_large,ii); -U_out(1:dim_n,ii)]);
%         end
%         kk                         = mm_real;
%         
%         for ii = 1:mm_cmp
%             [ U_in(:,kk+1:kk+2), U_out(:,kk+1:kk+2) ] = J_orthogonal( [real(ev_inside_c(:,idx_cmp_ew(ii))) ...
%                 imag(ev_inside_c(:,idx_cmp_ew(ii)))], [real(ev_all_outside_c(:,idx_cmp_ew(ii))) ...
%                 imag(ev_all_outside_c(:,idx_cmp_ew(ii)))], dim_n );
%             kk                                        = kk + 2;
%         end
% %         U_in(:,mm_real+1:2:mm   ) = real(ev_inside_c(:,idx_cmp_ew));
% %         U_in(:,mm_real+2:2:mm   ) = imag(ev_inside_c(:,idx_cmp_ew));
% %         
% %         U_out(:,        1:mm_real) = real(ev_all_outside_c(:,idx_real_ew));
% %         U_out(:,mm_real+1:2:mm   ) = real(ev_all_outside_c(:,idx_cmp_ew));
% %         U_out(:,mm_real+2:2:mm   ) = imag(ev_all_outside_c(:,idx_cmp_ew));
% % %         U_in        = [ real(ev_inside_c(:,idx_real_ew)), real(ev_inside_c(:,idx_cmp_ew)), imag(ev_inside_c(:,idx_cmp_ew)) ];
% % %         U_out       = [ real(ev_all_outside_c(:,idx_real_ew)), real(ev_all_outside_c(:,idx_cmp_ew)), imag(ev_all_outside_c(:,idx_cmp_ew)) ];
% %         
% %         for ii = 1:size(U_in,2)
% %             U_in(:,ii) = U_in(:,ii) / norm(U_in(:,ii));
% %         end
% %              
% %         for ii = 1:size(U_out,2)
% %             scale       = U_in(:,ii)' * [ U_out(dim_n+1:dim_large,ii); -U_out(1:dim_n,ii) ];
% %             U_out(:,ii) = U_out(:,ii) / scale;
% %         end
%         
% %         rsdl_test = zeros(mm_cmp,1);
% %         for ii = 1:mm_cmp
% %             tmp      = ev_all_outside_c(:,idx_cmp_ew(ii)) / norm(ev_all_outside_c(:,idx_cmp_ew(ii)));
% %             lambda   = 1/ew_inside_c(idx_cmp_ew(ii));
% %             rsdl_vec = mtx_M * tmp - lambda * mtx_L * tmp;
% %             rsdl_test(ii,1) = norm(rsdl_vec) / abs(lambda);
% %         end
% %         
% %         Err = U_in' * [ U_in(dim_n+1:dim_large,:); - U_in(1:dim_n,:) ];
% % 
% %         nrm_1 = zeros(size(Err,2),1);
% %         for ii = 1:size(Err,2)
% %             nrm_1(ii,1) = norm(Err(:,ii),inf);
% %         end
% %         
% %         Err2 = U_out' * [ U_out(dim_n+1:dim_large,:); - U_out(1:dim_n,:) ];
% % 
% %         nrm_2 = zeros(size(Err2,2),1);
% %         for ii = 1:size(Err2,2)
% %             nrm_2(ii,1) = norm(Err2(:,ii),inf);
% %         end
% %         
% %         Err3 = U_out' * [ U_in(dim_n+1:dim_large,:); - U_in(1:dim_n,:) ] + eye(size(U_out,2));
% %         nrm_3 = zeros(size(Err3,2),1);
% %         for ii = 1:size(Err3,2)
% %             nrm_3(ii,1) = norm(Err3(:,ii),inf);
% %         end
%         
%         U_on_u = [real(ev_inside_c(:,idx_unit_ew)), imag(ev_inside_c(:,idx_unit_ew))];
%         for ii = 1:mm_unit_ew
%             U_on_u(:,ii)            = U_on_u(:,ii) / norm(U_on_u(:,ii));
%             U_on_u(:,ii+mm_unit_ew) = U_on_u(:,ii+mm_unit_ew) / (U_on_u(:,ii)'*[U_on_u(dim_n+1:dim_large,ii+mm_unit_ew); -U_on_u(1:dim_n,ii+mm_unit_ew)]);
%         end
%         
% %         Err4 = U_on_u' * [ U_on_u(dim_n+1:dim_large,:); -U_on_u(1:dim_n,:)] - [zeros(mm_unit_ew), eye(mm_unit_ew); -eye(mm_unit_ew), zeros(mm_unit_ew)];
% 
%         U_eig_space = [ U_in U_on_u(:,1:mm_unit_ew) U_out  U_on_u(:,1+mm_unit_ew:2*mm_unit_ew) ];
%         X_init      = randn(dim_large, 2*(dim_n-no_cong_ew));
%         mtx_Y       = U_eig_space' * [ X_init(dim_n+1:dim_large,:); -X_init(1:dim_n,:) ];
%         X_new       = X_init - U_eig_space * [ -mtx_Y(no_cong_ew+1:2*no_cong_ew,:); mtx_Y(1:no_cong_ew,:) ];
%         
% %         Err         = U_eig_space' * [ X_new(dim_n+1:dim_large,:); -X_new(1:dim_n,:) ];
% %         fprintf('Err_J_orth = %11.4e \n', norm(Err, inf));
%         
%         RR_A = X_new' * (mtx_M * X_new);
%         RR_B = X_new' * (mtx_L * X_new);
%         
%         [RR_EV,RR_ew]  = eig(RR_A, RR_B, 'qz');
%         RR_ew          = diag(RR_ew);
%         idx_ew         = find(abs(RR_ew) <= 1+1.0e-9 & imag(RR_ew) > -1.0e-7);
%         RR_ew          = RR_ew(idx_ew);
%         RR_EV          = RR_EV(:,idx_ew);
%         
%         angle_ew = angle(RR_ew);
%         idx_tmp  = find( abs(angle_ew) <= 1.0e-10 );
%         if ( ~isempty(idx_tmp) )
%             angle_ew(idx_tmp) = 0;
%         end
%         angle_ew = mod(angle_ew,2*pi);
%         [~,idx_sort] = sort(angle_ew);
%         angle_ew     = angle_ew(idx_sort);
%         RR_ew        = RR_ew(idx_sort);
%         RR_EV        = RR_EV(idx_sort);
%         
%         kk           = 1;
%         while ( kk <= 20 && ~isempty(RR_ew) )
%             idx_tmp      = find(abs(angle_ew(1) - angle_ew(2:end)) < 15/180*pi);
%             if ( ~isempty(idx_tmp) )
%                 target_ew     = [RR_ew(1); RR_ew(idx_tmp)];
%                 idx_remainder = setdiff(1:length(RR_ew), [1; idx_tmp+1]);
%             else
%                 target_ew     = RR_ew(1);
%                 idx_remainder = 2:length(RR_ew);
%             end
%         
%             shift_S_Sinv  = sum(target_ew + 1./target_ew) / length(target_ew);
%             tmp_sqrt_root = sqrt(shift_S_Sinv^2 - 4);
%             tmp           = [ (shift_S_Sinv+tmp_sqrt_root)/2; (shift_S_Sinv-tmp_sqrt_root)/2 ];
%             [~,idx_min]   = min(abs(tmp - RR_ew(1)));
%             new_shift     = tmp(idx_min);
%         
%             %idx_remainder = setdiff(1:length(RR_ew), [1; idx_tmp+1]);
%             RR_ew         = RR_ew(idx_remainder);
%             RR_EV         = RR_EV(:,idx_remainder);
%             kk            = kk + 1;
%         end
%         
%     else
%         fprintf('Error !! no_cong_ew - dim_n = %2.0f \n', no_cong_ew - dim_n)
%     end
% end

     no_circle = no_circle - 1;
     
     switch Eng_eps 
        case 0
            save tmp_dat_eps0 rad rad_old ew_unit_c ev_on_c ew_inside_c no_outside_ev ...
                ev_inside_c ref_ev_outside_c no_circle ref_ew_inside_c ref_ev_inside_c ...
                max_wanted_ew ew_inner rad_circle radius_case weighting_choosing ev_all_outside_c
            
        case 0.5
            save tmp_dat_eps0p5 rad rad_old ew_unit_c ev_on_c ew_inside_c no_outside_ev ...
                ev_inside_c ref_ev_outside_c no_circle ref_ew_inside_c ref_ev_inside_c ...
                max_wanted_ew ew_inner rad_circle radius_case weighting_choosing ev_all_outside_c
            
        case 1
            save tmp_dat_eps1p0 rad rad_old  ew_unit_c ev_on_c ew_inside_c no_outside_ev ...
                ev_inside_c ref_ev_outside_c no_circle ref_ew_inside_c ref_ev_inside_c ...
                max_wanted_ew ew_inner rad_circle radius_case weighting_choosing ev_all_outside_c
            
        case 1.5
            save tmp_dat_eps1p5 rad rad_old ew_unit_c ev_on_c ew_inside_c no_outside_ev ...
                ev_inside_c ref_ev_outside_c no_circle ref_ew_inside_c ref_ev_inside_c ...
                max_wanted_ew ew_inner rad_circle radius_case weighting_choosing ev_all_outside_c
            
        case 2
            save tmp_dat_eps2p0 rad rad_old ew_unit_c ev_on_c ew_inside_c no_outside_ev...
                ev_inside_c ref_ev_outside_c no_circle ref_ew_inside_c ref_ev_inside_c ...
                max_wanted_ew ew_inner rad_circle radius_case weighting_choosing ev_all_outside_c
            
        case 2.5
            save tmp_dat_eps2p5 rad rad_old ew_unit_c ev_on_c ew_inside_c no_outside_ev ...
                ev_inside_c ref_ev_outside_c no_circle ref_ew_inside_c ref_ev_inside_c ...
                max_wanted_ew ew_inner rad_circle radius_case weighting_choosing ev_all_outside_c
            
        case 3
            save tmp_dat_eps3p0 rad rad_old ew_unit_c ev_on_c ew_inside_c no_outside_ev ...
                ev_inside_c ref_ev_outside_c no_circle ref_ew_inside_c ref_ev_inside_c ...
                max_wanted_ew ew_inner rad_circle radius_case weighting_choosing ev_all_outside_c
            
        case 3.5
            save tmp_dat_eps3p5 rad rad_old ew_unit_c ev_on_c ew_inside_c no_outside_ev ...
                ev_inside_c ref_ev_outside_c no_circle ref_ew_inside_c ref_ev_inside_c ...
                max_wanted_ew ew_inner rad_circle radius_case weighting_choosing ev_all_outside_c
            
        case 4
            save tmp_dat_eps4p0 rad rad_old ew_unit_c ev_on_c ew_inside_c no_outside_ev ...
                ev_inside_c ref_ev_outside_c no_circle ref_ew_inside_c ref_ev_inside_c ...
                max_wanted_ew ew_inner rad_circle radius_case weighting_choosing ev_all_outside_c
            
        case 4.5
            save tmp_dat_eps4p5 rad rad_old ew_unit_c ev_on_c ew_inside_c no_outside_ev ...
                ev_inside_c ref_ev_outside_c no_circle ref_ew_inside_c ref_ev_inside_c ...
                max_wanted_ew ew_inner rad_circle radius_case weighting_choosing ev_all_outside_c
            
        case 5
            save tmp_dat_eps5p0 rad rad_old ew_unit_c ev_on_c ew_inside_c no_outside_ev ...
                ev_inside_c ref_ev_outside_c no_circle ref_ew_inside_c ref_ev_inside_c ...
                max_wanted_ew ew_inner rad_circle radius_case weighting_choosing ev_all_outside_c
    
        case 5.5
            save tmp_dat_eps5p5 rad rad_old ew_unit_c ev_on_c ew_inside_c no_outside_ev ...
                ev_inside_c ref_ev_outside_c no_circle ref_ew_inside_c ref_ev_inside_c ...
                max_wanted_ew ew_inner rad_circle radius_case weighting_choosing ev_all_outside_c
            
        case 6
            save tmp_dat_eps6p0 rad rad_old ew_unit_c ev_on_c ew_inside_c no_outside_ev ...
                ev_inside_c ref_ev_outside_c no_circle ref_ew_inside_c ref_ev_inside_c ...
                max_wanted_ew ew_inner rad_circle radius_case weighting_choosing ev_all_outside_c
            
        case 7
            save tmp_dat_eps7p0 rad rad_old ew_unit_c ev_on_c ew_inside_c no_outside_ev ...
                ev_inside_c ref_ev_outside_c no_circle ref_ew_inside_c ref_ev_inside_c ...
                max_wanted_ew ew_inner rad_circle radius_case weighting_choosing ev_all_outside_c
            
         case -2
             save tmp_dat_epsn2p0 rad rad_old ew_unit_c ev_on_c ew_inside_c no_outside_ev ...
                ev_inside_c ref_ev_outside_c no_circle ref_ew_inside_c ref_ev_inside_c ...
                max_wanted_ew ew_inner rad_circle radius_case weighting_choosing ev_all_outside_c
            
         case -5.1
             save tmp_dat_epsn5p1 rad rad_old ew_unit_c ev_on_c ew_inside_c no_outside_ev ...
                ev_inside_c ref_ev_outside_c no_circle ref_ew_inside_c ref_ev_inside_c ...
                max_wanted_ew ew_inner rad_circle radius_case weighting_choosing ev_all_outside_c
            
         case -6.5
             save tmp_dat_epsn6p5 rad rad_old ew_unit_c ev_on_c ew_inside_c no_outside_ev ...
                ev_inside_c ref_ev_outside_c no_circle ref_ew_inside_c ref_ev_inside_c ...
                max_wanted_ew ew_inner rad_circle radius_case weighting_choosing ev_all_outside_c
            
         otherwise
             save tmp_dat rad rad_old ew_unit_c ev_on_c ew_inside_c no_outside_ev ...
                ev_inside_c ref_ev_outside_c no_circle ref_ew_inside_c ref_ev_inside_c ...
                max_wanted_ew ew_inner rad_circle radius_case weighting_choosing ev_all_outside_c
            
     end

end

if ( strcmp(plot_fig, 'yes') )
    figure(2)
    hold off
    close(vidObj);

    figure(3)
    hold off 

    figure(1)
    hold off
end

end
%diary off
