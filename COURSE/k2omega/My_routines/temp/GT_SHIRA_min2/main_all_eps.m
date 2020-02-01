clear
format long e
%clc

% Eng_eps = 5.5; %0.5; %1.0; %0; %1.5; %7; %5.0; %4.5;

load ../0804/h00L
mtx_H0_large = A;
load ../0804/h01L
mtx_H1_large = A;

load ../0804/s00L
mtx_S0_large = A;
load ../0804/s01L
mtx_S1_large = A;

dim_large = size(mtx_H0_large,1);
dim_n     = dim_large / 2;
mtx_H0    = mtx_H0_large(1:dim_n, 1:dim_n);
mtx_H1    = mtx_H1_large(1:dim_n,dim_n+1:2*dim_n);
mtx_S0    = mtx_S0_large(1:dim_n, 1:dim_n);
mtx_S1    = mtx_S1_large(1:dim_n,dim_n+1:2*dim_n);

ij                  = 1;
Eps_interval        = -20:0.1:10;
no_unit_ew          = zeros(length(Eps_interval),1);
iter_all            = zeros(length(Eps_interval),1);
CPUTime_all         = zeros(length(Eps_interval),1);
ratio_deuplicate_ew = zeros(length(Eps_interval),1);

gamma            = 1; %exp(pi/12*1i); %[1; exp(1i*(theta(2:end-1)'+ 1.0e-3)); -1]; %1;%exp(1i*(theta(2:end-1)'+ 1.0e-3));
leng_r           = length(gamma);
CPUTime_LU       = zeros(leng_r,1);
CPUtime_JLanczos = zeros(leng_r,1);
CPUtime_gen_M    = zeros(leng_r,1);
CPUtime_eig      = zeros(leng_r,1);
NO_ew            = zeros(leng_r,1);
Es_ew_symplectic = zeros(4*dim_n,leng_r);
ew_Hamiltonian   = zeros(2*dim_n,leng_r);

plot_flag = 'no';

theta = 0:0.001:2*pi;
x     = cos(theta)+1i*sin(theta);

for Eng_eps = Eps_interval
    fprintf('eps = %11.4e \n', Eng_eps);
    mtx_A   = Eng_eps * mtx_S1 - mtx_H1;
    mtx_Q   = Eng_eps * mtx_S0 - mtx_H0;

    mtx_NME.mtx_A = Eng_eps * mtx_S1 - mtx_H1;
    mtx_NME.mtx_Q = Eng_eps * mtx_S0 - mtx_H0;
    
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
     
    mtx_L = sparse(i_idx_L, j_idx_L, val_L);

    figure(1) 
    plot(x,'r-')
    axis([-1.1, 1.1, -1.1, 1.1]);
    hold on

    figure(2) 
    plot(x,'r-')
    axis([-1.1, 1.1, -1.1, 1.1]);
    title(['eps = ' num2str(Eng_eps)]);
    hold on

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
        res = zeros(2*dim_n,1);
        for ii = 1:2*dim_n
            res(ii,1) = norm(Err(:,ii));
        end
    
        save EW_data_n5p1 ew res
    
    else
    
        switch Eng_eps 
            case 0
                load ../0804/EW_data_0
            
            case 0.5
                load ../0804/EW_data_0p5
            
            case 1
                load ../0804/EW_data_1p0
            
            case 1.5
                load ../0804/EW_data
            
            case 2
                load ../0804/EW_data_2p0
            
            case 2.5
                load ../0804/EW_data_2p5
            
            case 3
                load ../0804/EW_data_3p0
            
            case 3.5
                load ../0804/EW_data_3p5
            
            case 4
                load ../0804/EW_data_4p0
            
            case 4.5
                load ../0804/EW_data_4p5
            
            case 5
                load ../0804/EW_data_5p0
    
            case 5.5
                load ../0804/EW_data_55
            
            case 6
                load ../0804/EW_data_60
            
            case 7
                load ../0804/EW_data_7p0
                
            case -2
                load ../v1_0919/EW_data_n2p0 
            
            case -5.1
                load ../v2_0924/EW_data_n5p1 
            
            case -6.5
                load ../v1_0919/EW_data_n6p5 
            
        end
    end

    plot_eig = 'no';

    if ( strcmp(plot_eig, 'yes') )
        figure(1)
        plot(ew,'bo','LineWidth',1); 
        hold on

        figure(2)
        plot(ew,'bo','LineWidth',1); 
        hold on
    end

    deflation = 'yes'; %'yes';
    switch deflation 
        case 'yes'
            deflation_force = 'yes';
            
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
            vidObj           = [];
            
            tic;
            
            [ ew_unit_c, ev_on_c, ew_inside_c, ev_inside_c, no_deuplicate_ew, ...
                rad, ev_outside_c, ref_ew_inside_c, ref_ev_inside_c, max_wanted_ew, flag_conv, ...
                ew_inner, radius_case, weighting_choosing, iter_all(ij,1) ] = compute_1st_circle_ew( dim_n, mtx_NME, ...
                mtx_M, mtx_L, rad, rad_old, max_wanted_ew, radius_case, vidObj, weighting_choosing, ...
                ew_unit_c, ev_on_c, ew_inside_c, ev_inside_c, ev_outside_c, ref_ew_inside_c, ref_ev_inside_c, plot_flag );
            
            CPUTime_all(ij,1)         = toc;
            ratio_deuplicate_ew(ij,1) = sum(no_deuplicate_ew(:,1)) / sum(no_deuplicate_ew(:,2));
            
%             tic;
%             [ ew_unit_c, ev_on_c, ew_inside_c, ev_inside_c, ev_outside_c, no_deuplicate_ew, iter_all(ij,1), ...
%                 rad, exceptional, rad_check, ew_inner ] = compute_ew_near_unit_deflate_1010( dim_n, mtx_NME, deflation_force );
%             CPUTime_all(ij,1)         = toc;
%             ratio_deuplicate_ew(ij,1) = sum(no_deuplicate_ew(:,1)) / sum(no_deuplicate_ew(:,2));

        case 'no'
            deflation_force = 'no';
            
            tic;
            [ ew_unit_c, ev_on_c, ew_inside_c, ev_inside_c, ev_outside_c, no_deuplicate_ew, iter_all(ij,1), ...
                rad, exceptional, rad_check, ew_inner ] = compute_ew_near_unit_deflate_1010( dim_n, mtx_NME, deflation_force );
            CPUTime_all(ij,1)         = toc;
            ratio_deuplicate_ew(ij,1) = sum(no_deuplicate_ew(:,1)) / sum(no_deuplicate_ew(:,2));
            
    end
   
    no_unit_ew(ij) = 2 * size(ew_unit_c,1);
    ij             = ij + 1;
   
    figure(1)
    hold off
   
    figure(2)
    hold off
   
    figure(3)
    plot(Eps_interval(1:ij-1), no_unit_ew(1:ij-1),'bo','LineWidth',2);
    
    switch deflation
        case 'yes'
            save all_eps_unit_ew_defl_neg20p0_10p0 Eps_interval no_unit_ew CPUTime_all iter_all ratio_deuplicate_ew
    
        case 'no'
            save all_eps_unit_ew_without_defl Eps_interval no_unit_ew CPUTime_all iter_all ratio_deuplicate_ew
    end
end

figure(1)
hold off
