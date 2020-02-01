function vec_z = FAME_Matrix_Vector_Production_IFFT_Triple_General(vec_x,Nx,Ny,Nz,N,D_kx,F_ky,F_kz)
    
    vec_z = zeros(3*N,1);
    for i = 1:3
        vec_x_temp = reshape(vec_x((i-1)*N+1:i*N), Nz, Nx*Ny);
%         vec_x_temp = conj(D_kz).*(ifft(vec_x_temp));
        vec_x_temp = F_kz.*(ifft(vec_x_temp));

        vec_y_temp = zeros(Ny, Nx*Nz);

        for ii = 1:Nx
            vec_y_temp(:,(ii-1)*Nz+1:ii*Nz) = vec_x_temp(:,(ii-1)*Ny+1:ii*Ny).';
        end

        vec_y_temp = ifft(vec_y_temp);

        idx_col = zeros(Nx*Nz,1);
        for ii = 1:Nz
            idx_col((ii-1)*Nx+1:ii*Nx,1) = (ii:Nz:(Nx-1)*Nz+ii)';
        end

%         vec_y_temp = conj(D_ky)*vec_y_temp(:,idx_col);
        vec_y_temp = F_ky.*vec_y_temp(:,idx_col);

        vec_z_temp = zeros(Nx, Ny*Nz);
        for ii = 1:Nz
            vec_z_temp(:,(ii-1)*Ny+1:ii*Ny) = vec_y_temp(:,(ii-1)*Nx+1:ii*Nx).';
        end

        % mtx_Q_x_tilde = ifft(mtx_Q_x_tilde);

%         vec_z((i-1)*N+1:i*N) = sqrt(Nx * Ny * Nz) * reshape(sparse(1:Nx,1:Nx,conj(D_kx))*ifft(vec_z_temp),Nx*Ny*Nz,1);
        vec_z((i-1)*N+1:i*N) = sqrt(Nx * Ny * Nz) * reshape(sparse(1:Nx,1:Nx,D_kx)*ifft(vec_z_temp),Nx*Ny*Nz,1);
    end
end
