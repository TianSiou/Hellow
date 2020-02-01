function [ freq, Ele_field, cpu_time, LS_iter, LS_cpu_time ] = FAME_Fast_Algorithms_Simple_Anisotropic_fem(grid_num, wave_vec, B, Lambdas, eigen_wanted )
    global inner_iter inner_count inner_cpu_time
    inner_count = 0;  
    
    Nx = grid_num(1); Ny = grid_num(2); Nz = grid_num(3);
    N  = Nx*Ny*Nz;
    if norm(wave_vec) == 0
        Nd = 2*N - 2;
    else 
        Nd = 2*N;
    end
    Lambdas.D_k   = kron(Lambdas.D_kz,kron(Lambdas.D_ky,Lambdas.D_kx));
    Lambdas.D_ks  = conj(Lambdas.D_k);
    
    opt.isreal = 0;
    opt.issym  = 1;
    opt.tol    = 1e-12;
    LS_tol     = 1e-12;
    
    time_eig = tic;
    [ev,ew] = eigs( @(x) FAME_Matrix_Vector_Production_Anisotropic_invAr_Simple_fem( x, B, Nx, Ny, Nz, N, Lambdas.Sigma_r, Lambdas.Pi_Qr, Lambdas.Pi_Qrs, Lambdas.D_k, Lambdas.D_ks, LS_tol),...
                         Nd, eigen_wanted, 'lm', opt );
    cpu_time = toc(time_eig);
    
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
        
        Output_eigvec = FAME_Matrix_Vector_Production_invB_Anisotropic_fem(vec_y, B);
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
    [Output_eigval,idx] = sort(Output_eigval,'descend');
    Output_eigvec_mat = Output_eigvec_mat(:,idx);
end