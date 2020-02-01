function [ freq, Ele_field, cpu_time, LS_iter, LS_cpu_time ] = FAME_Fast_Algorithms_Simple_Isotropic_gpu( grid_num, wave_vec, B, Lambdas, eigen_wanted )
    global inner_iter inner_count inner_cpu_time
    inner_count = 0;
    
    Nx = grid_num(1); Ny = grid_num(2); Nz = grid_num(3);
    
    N = Nx*Ny*Nz;
    if norm(wave_vec) == 0
        Nd = 2*N - 2;
    else 
        Nd = 2*N;
    end
    fprintf('start h2d...')
    tic
    [dB, dLambdas] = Host2Device(B,Lambdas);
    % test
    vec_x  = rand(Nd,1);
    dvec_x = gpuArray(vec_x);
    time_h2d = toc
    for i = 1:100
        tic
        dvec_y = FAME_Matrix_Vector_Production_Isotropic_invAr_Simple( dvec_x, dB, Nx, Ny, Nz, N, dLambdas.Sigma_r, dLambdas.Pi_Qr, dLambdas.Pi_Qrs, dLambdas.D_k, dLambdas.D_ks, 1e-12 );
%         dvec_y = FAME_Matrix_Vector_Production_Isotropic_Ar_Simple(dvec_x, dB, Nx, Ny, Nz, N, dLambdas.Pi_Qr, dLambdas.Pi_Qrs, dLambdas.D_k, dLambdas.D_ks);
%         dvec_y = FAME_Matrix_Vector_Production_Isotropic_Ar_Simple(dvec_x, dB, Nx, Ny, Nz, N, dLambdas.Pi_Qr, dLambdas.Pi_Qrs, dLambdas.D_k, dLambdas.D_ks);
        time_d = toc
        tic
        vec_y = FAME_Matrix_Vector_Production_Isotropic_invAr_Simple( vec_x, B, Nx, Ny, Nz, N, Lambdas.Sigma_r, Lambdas.Pi_Qr, Lambdas.Pi_Qrs, Lambdas.D_k, Lambdas.D_ks, 1e-12 );
%         vec_y = FAME_Matrix_Vector_Production_Isotropic_Ar_Simple(vec_x, B, Nx, Ny, Nz, N, Lambdas.Pi_Qr, Lambdas.Pi_Qrs, Lambdas.D_k, Lambdas.D_ks);
%         vec_y = FAME_Matrix_Vector_Production_Isotropic_Ar_Simple(vec_x, B, Nx, Ny, Nz, N, Lambdas.Pi_Qr, Lambdas.Pi_Qrs, Lambdas.D_k, Lambdas.D_ks);
        time_h = toc
    end
    %=====================
    
    
    opt.isreal = 0;
    opt.issym  = 1;
    
    LS_tol = 1e-12;
    
    cpu_time_start = tic;
    [ev,ew] = eigs( @(x) FAME_Matrix_Vector_Production_Isotropic_invAr_Simple( x, B, Nx, Ny, Nz, N, Lambdas.Sigma_r, Lambdas.Pi_Qr, Lambdas.Pi_Qrs, Lambdas.D_k, Lambdas.D_ks, LS_tol),...
                         Nd, eigen_wanted, 'lm', opt );
    cpu_time = toc(cpu_time_start);
    
    LS_iter     = inner_iter;
    LS_cpu_time = inner_cpu_time;
    
    clear('inner_iter','inner_cpu_time');
    
    ew = 1./diag(ew);
    [ Ele_field, freq ] = Eigen_Restoration_Simple_Isotropic( ev, ew, Nx, Ny, Nz, N, wave_vec, B, Lambdas );
end

function [ Output_eigvec_mat, Output_eigval ] = Eigen_Restoration_Simple_Isotropic( Input_eigvec_mat, Input_eigval, Nx, Ny, Nz, N, wave_vec, B, Lambdas )
%% Start to restore the eigenvectors
    for column_idx = 1 : size(Input_eigvec_mat,2)
        vec_y = Lambdas.Sigma_r.*Input_eigvec_mat(:,column_idx);
        vec_y = FAME_Matrix_Vector_Production_Qr_Simple(vec_y, Nx, Ny, Nz, N, Lambdas.Pi_Qr, Lambdas.Pi_Qrs, Lambdas.D_k, Lambdas.D_ks, 'normal');

%         Output_eigvec = invB_eps.*vec_y;     
        Output_eigvec = FAME_Matrix_Vector_Production_invB_Isotropic(vec_y, B);
% Normalizing the eigenvector     
        Output_eigvec = Output_eigvec / norm(Output_eigvec);           
        Output_eigvec_mat(:,column_idx) = Output_eigvec;        
    end
    Output_eigval = Input_eigval;
    if  norm(wave_vec) == 0  
        Output_eigvec_mat = [ ones(3*N,2), Output_eigvec_mat ];
        Output_eigval     = [0;0;Output_eigval];
        Output_eigvec_mat = Output_eigvec_mat(:,1:end-2);
        Output_eigval     = Output_eigval(1:end-2);
    end
    Output_eigval = sqrt(Output_eigval);
    [Output_eigval,idx] = sort(real(Output_eigval));
    Output_eigvec_mat   = Output_eigvec_mat(:,idx);
end

function [dB, dLambdas] = Host2Device(B,Lambdas)
    dB.invB_eps   = gpuArray(B.invB_eps);
    dLambdas.D_k  = gpuArray(Lambdas.D_k);
    dLambdas.D_ks = gpuArray(Lambdas.D_ks);
    dLambdas.Pi_Qr  = gpuArray(Lambdas.Pi_Qr);
    dLambdas.Pi_Qrs = gpuArray(Lambdas.Pi_Qrs);
    dLambdas.Sigma_r = gpuArray(Lambdas.Sigma_r);
end