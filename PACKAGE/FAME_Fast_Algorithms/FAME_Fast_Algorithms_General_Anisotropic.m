function [ freq, Ele_field, cpu_time, LS_iter, LS_cpu_time ] = FAME_Fast_Algorithms_General_Anisotropic(grid_num, wave_vec, B, Lambdas, eigen_wanted)
    global inner_iter inner_count inner_cpu_time
    inner_count = 0;    

    Nx = grid_num(1); Ny = grid_num(2); Nz = grid_num(3);
    N  = Nx*Ny*Nz;
    if norm(wave_vec) == 0
        Nd = 2*N - 2;
    else 
        Nd = 2*N;
    end
    
    opt.isreal = 0;
    opt.issym  = 1;
    LS_tol = 1e-12;
    
    cpu_time_start = tic;
    [ev,ew] = eigs( @(x) FAME_Matrix_Vector_Production_Anisotropic_invAr_General( x, B, Nx, Ny, Nz, N,...
                         Lambdas.Sigma_r, Lambdas.Pi_Qr, Lambdas.Pi_Qrs, Lambdas.D_kx, Lambdas.D_ky, Lambdas.D_kz, Lambdas.F_ky, Lambdas.F_kz, LS_tol),...
                         Nd, eigen_wanted, 'lm', opt );
    cpu_time = toc(cpu_time_start);
    
    LS_iter     = inner_iter;
    LS_cpu_time = inner_cpu_time;
    
    clear('inner_iter','inner_cpu_time');
    
    ew = 1./diag(ew);
    [ Ele_field, freq ] = Eigen_Restoration_General_Isotropic( ev, ew, Nx, Ny, Nz, N, wave_vec, B, Lambdas );
end

function [ Output_eigvec_mat, Output_eigval ] = Eigen_Restoration_General_Isotropic( Input_eigvec_mat, Input_eigval, Nx, Ny, Nz, N, wave_vec, B, Lambdas )
%% Start to restore the eigenvectors
    for column_idx = 1 : size(Input_eigvec_mat,2)
        vec_y = Lambdas.Sigma_r.*Input_eigvec_mat(:,column_idx);
        vec_y = FAME_Matrix_Vector_Production_Qr_General(vec_y, Nx, Ny, Nz, N, Lambdas.Pi_Qr, Lambdas.Pi_Qrs, Lambdas.D_kx, Lambdas.D_ky, Lambdas.D_kz, Lambdas.F_ky, Lambdas.F_kz, 'normal');
        
%         Output_eigvec = [ invB_eps.d11.*vec_y(1:N) + invB_eps.d12.*vec_y(N+1:2*N) + invB_eps.d13.*vec_y(2*N+1:3*N);
%                           invB_eps.d21.*vec_y(1:N) + invB_eps.d22.*vec_y(N+1:2*N) + invB_eps.d23.*vec_y(2*N+1:3*N);
%                           invB_eps.d31.*vec_y(1:N) + invB_eps.d32.*vec_y(N+1:2*N) + invB_eps.d33.*vec_y(2*N+1:3*N) ];
        Output_eigvec = FAME_Matrix_Vector_Production_invB_Anisotropic(vec_y, B, N);
%         Output_eigvec = invB_eps*vec_y;       
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
end