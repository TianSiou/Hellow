function [ freq, Ele_field, Mag_field, cpu_time, LS_iter, LS_cpu_time ] = FAME_Fast_Algorithms_Simple_Biisotropic(grid_num, wave_vec, B, Lambdas, eigen_wanted)
    global inner_iter inner_count inner_cpu_time
    inner_count = 0;

    Nx = grid_num(1); Ny = grid_num(2); Nz = grid_num(3);
    
    N = Nx*Ny*Nz;
    if norm(wave_vec) == 0
        Nd = 4*N - 4;
    else 
        Nd = 4*N;
    end
    
%     Lambdas.D_k              = kron(Lambdas.D_kz,kron(Lambdas.D_ky,Lambdas.D_kx));
%     Lambdas.D_ks             = conj(Lambdas.D_k);
    B.B_zeta_s               = conj(B.B_zeta);
    B.invPhi                 = 1./( B.B_eps - (B.B_xi./B.B_mu).*B.B_zeta );
    Lambdas.inv_Sigma_r      = 1./ Lambdas.Sigma_r;
    
    eigs_opt.isreal = 0;
    eigs_opt.issym  = 0;
    
    LS_tol = 1e-12;
    
    cpu_time_start = tic;
    [ ev, ew ]     = eigs( @(x) FAME_Matrix_Vector_Production_Biisotropic_Simple_Posdef( x, B, Nx, Ny, Nz, N, Lambdas.inv_Sigma_r, Lambdas.Pi_Qr, Lambdas.Pi_Pr, Lambdas.Pi_Qrs, Lambdas.Pi_Prs, Lambdas.D_k, Lambdas.D_ks, LS_tol),...
                           Nd, eigen_wanted, 'lr', eigs_opt );
    cpu_time       = toc(cpu_time_start);

    LS_iter        = inner_iter;
    LS_cpu_time    = inner_cpu_time;
    
    clear('inner_iter','inner_cpu_time');
    
    ew = 1./diag(ew);
    [ EleMag_field, freq ] = Eigen_Restoration_Simple_Biisotropic( ev, ew, Nx, Ny, Nz, N, wave_vec, B, Lambdas);
    Ele_field = EleMag_field(1:end/2,:);
    Mag_field = EleMag_field(end/2+1:end,:);
    
end

function [ Output_eigvec_mat, Output_eigval ] = Eigen_Restoration_Simple_Biisotropic( Input_eigvec_mat, Input_eigval, Nx, Ny, Nz, N, wave_vec, B, Lambdas)
idx = find( real(Input_eigval)> -1e-6  );
Input_eigval = Input_eigval(idx);
Input_eigvec_mat = Input_eigvec_mat(:,idx);
%% Start to restore the eigenvectors
    for column_idx = 1 : size(Input_eigvec_mat,2)
        vec_x     = Input_eigvec_mat(:,column_idx);
        vec_x_ele = vec_x(1:end/2);
        vec_x_mag = vec_x(end/2+1:end);   
        
        vec_y_ele = FAME_Matrix_Vector_Production_Pr_Simple(vec_x_ele, Nx, Ny, Nz, N, Lambdas.Pi_Pr, Lambdas.Pi_Prs, Lambdas.D_k, Lambdas.D_ks, 'normal');
        vec_y_mag = FAME_Matrix_Vector_Production_Qr_Simple(vec_x_mag, Nx, Ny, Nz, N, Lambdas.Pi_Qr, Lambdas.Pi_Qrs, Lambdas.D_k, Lambdas.D_ks, 'normal');
        vec_y     = [vec_y_ele;vec_y_mag];
        
%         Output_eigvec = B\vec_y;    
        Output_eigvec = FAME_Matrix_Vector_Production_invB_Biisotropic(vec_y, B);
% Normalizing the eigenvector     
        Output_eigvec = Output_eigvec / norm(Output_eigvec);           
        Output_eigvec_mat(:,column_idx) = Output_eigvec;        
    end
    Output_eigval = Input_eigval;
    if  norm(wave_vec) == 0  
        [~,idx] = sort(abs(Output_eigval),'descend');
        Output_eigval(idx(1:2)) = 0;
        Output_eigvec_mat(:,idx(1:2)) = ones(6*N,2);
    end
    Output_eigval = real(Output_eigval);
    [Output_eigval,idx] = sort(Output_eigval);
    Output_eigvec_mat = Output_eigvec_mat(:,idx);
end