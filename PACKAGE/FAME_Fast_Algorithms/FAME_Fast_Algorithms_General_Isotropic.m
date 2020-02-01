function [ freq, Ele_field, cpu_time, LS_iter, LS_cpu_time ,prv] = FAME_Fast_Algorithms_General_Isotropic( grid_num, wave_vec, B, Lambdas, eigen_wanted )
    global inner_iter inner_count inner_cpu_time
    inner_count = 0;
    inner_iter = [];
    inner_cpu_time = [];
    shift = 0 * pi^2;
    
    
    Nx = grid_num(1); Ny = grid_num(2); Nz = grid_num(3);
    N = Nx*Ny*Nz;
    if norm(wave_vec) == 0
        Nd = 2*N - 2;
    else 
        Nd = 2*N;
    end

%     invB_eps = 1./B_eps;
    
    opt.isreal = 0;
    opt.issym  = 1;
    
    LS_tol = 1e-12;
    
    cpu_time_start = tic;
    
% LSopt.tol         = 1.0e-3; %1e-3
% LSopt.infoConvBeh = 'yes';    
%  stop_tolerance = 1.0e-12; %1.0e8 * eps / scale; 
%         target_type    = 'RGTR'; %'RGTC'; 
% 
%         no_restart     = 24; %35;
%         CorEq          = 'SIRA'; %'JD'; %'SIRA';
%         LSinfo.solver  = 'pcg';%'bicgstabl'; %'bicgstabl'; %'minres';
%         LSinfo.precond = 'no'; %'yes';
%         flag_LU = 0;
%         
%         if size(prv,1) == 0
%             initial_V   = randn(Nd,1) + 1i * randn(Nd,1); % ones(dim,1);
%             initial_V   = initial_V / norm(initial_V);
%         else
%             initial_V = prv;
%         end
%         
%         fid = 0;
%         target = 0*pi^2;
%      
%         [ew, ev] = GEP_AB_Herm_JDSIRA_Driver (@(x) FAME_Matrix_Vector_Production_Isotropic_Ar_General( x, B, Nx, Ny, Nz, N,...
%                           Lambdas.Pi_Qr, Lambdas.Pi_Qrs, Lambdas.D_kx, Lambdas.D_ky, Lambdas.D_kz, Lambdas.F_ky, Lambdas.F_kz),...
%                          @(x) (Lambdas.Sigma_r.^2) .\ x, Nd, no_restart, ...
%             eigen_wanted, stop_tolerance, initial_V, target, fid, target_type, CorEq, LSinfo); % @(x)solve_Minv_b( x, LU_precond ));             
%         [ew, i] = sort(ew);
%         ev      = ev(:,i); 
%         ev      = Lambdas.Sigma_r .\ ev;
%         prv     = ev(:, 1:5);
    if shift == 0
    [ev,ew] = eigs( @(x) FAME_Matrix_Vector_Production_Isotropic_invAr_General( x, B, Nx, Ny, Nz, N,...
                         Lambdas.Sigma_r, Lambdas.Pi_Qr, Lambdas.Pi_Qrs, Lambdas.D_kx, Lambdas.D_ky, Lambdas.D_kz, Lambdas.F_ky, Lambdas.F_kz, LS_tol),...
                         Nd, eigen_wanted, 'lm', opt );
    else
    [ev,ew] = eigs( @(x) FAME_Matrix_Vector_Production_Isotropic_shift_invAr_General( x, B, Nx, Ny, Nz, N,...
                         Lambdas.Sigma_r, Lambdas.Pi_Qr, Lambdas.Pi_Qrs, Lambdas.D_kx, Lambdas.D_ky, Lambdas.D_kz, Lambdas.F_ky, Lambdas.F_kz, LS_tol, shift),...
                         Nd, eigen_wanted, 'lm', opt );    
    end
    
    cpu_time = toc(cpu_time_start);
    
    LS_iter     = inner_iter;
    LS_cpu_time = inner_cpu_time;
    
    clear('inner_iter','inner_cpu_time');
   
    ew = diag(ew);
    ew = 1./ew + shift;
    
    [ Ele_field, freq ] = Eigen_Restoration_General_Isotropic( ev, ew, Nx, Ny, Nz, N, wave_vec, B, Lambdas, shift);
end

function [ Output_eigvec_mat, Output_eigval ] = Eigen_Restoration_General_Isotropic( Input_eigvec_mat, Input_eigval, Nx, Ny, Nz, N, wave_vec, B, Lambdas, shift)
%% Start to restore the eigenvectors
    for column_idx = 1 : size(Input_eigvec_mat,2)
        vec_y = Lambdas.Sigma_r.*Input_eigvec_mat(:,column_idx);
        vec_y = FAME_Matrix_Vector_Production_Qr_General(vec_y, Nx, Ny, Nz, N, Lambdas.Pi_Qr, Lambdas.Pi_Qrs, Lambdas.D_kx, Lambdas.D_ky, Lambdas.D_kz, Lambdas.F_ky, Lambdas.F_kz, 'normal');

%         Output_eigvec = invB_eps.*vec_y;   
        Output_eigvec = FAME_Matrix_Vector_Production_invB_Isotropic(vec_y, B);
% Normalizing the eigenvector     
        Output_eigvec = Output_eigvec / norm(Output_eigvec);           
        Output_eigvec_mat(:,column_idx) = Output_eigvec;        
    end
    Output_eigval = Input_eigval;
    if shift == 0
        if  norm(wave_vec) == 0  
            Output_eigvec_mat = [ ones(3*N,2), Output_eigvec_mat ];
            Output_eigval     = [0;0;Output_eigval];
            Output_eigvec_mat = Output_eigvec_mat(:,1:end-2);
            Output_eigval     = Output_eigval(1:end-2);
        end
    end
    Output_eigval = sqrt(Output_eigval);
    [Output_eigval,idx] = sort(real(Output_eigval));
    Output_eigvec_mat   = Output_eigvec_mat(:,idx);
end