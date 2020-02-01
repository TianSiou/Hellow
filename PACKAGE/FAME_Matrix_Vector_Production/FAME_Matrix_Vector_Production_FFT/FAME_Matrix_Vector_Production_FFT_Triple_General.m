function vec_z = FAME_Matrix_Vector_Production_FFT_Triple_General(vec_x,Nx,Ny,Nz,N,D_kx,D_ky,D_kz)
    vec_z = zeros(3*N,1);
    for i = 1:3
        vec_x_temp = fft(sparse(1:Nx,1:Nx,conj(D_kx)) * reshape(vec_x((i-1)*N+1:i*N), Nx, Nz*Ny));

        vec_y_temp = zeros(Ny, Nx*Nz);
        for ii = 1:Nx
            vec_y_temp(:,(ii-1)*Nz+1:ii*Nz) = sparse(1:Ny,1:Ny,conj(D_ky(:,ii))) * reshape(vec_x_temp(ii,:).', Ny, Nz);
        %     mtx_P_y_tilde(:,(ii-1)*Nz+1:ii*Nz) = sparse(1:Ny,1:Ny,conj(FFT_parameter.mtx_D_jx(:,ii))) * reshape(mtx_P_x_tilde(ii,:).', Ny, Nz);
        end

        vec_y_temp = fft(vec_y_temp);

        vec_z_temp = zeros(Nz, Nx*Ny);
        for ii = 1:Nx
            vec_z_temp(:,(ii-1)*Ny+1:ii*Ny) = conj(D_kz(:,1:Ny,ii)).*(vec_y_temp(:,(ii-1)*Nz+1:ii*Nz).');
%             mtx_P_z_tilde(:,(ii-1)*Ny+1:ii*Ny) = conj(FFT_parameter_bar.mtx_D_jell(:,1:Ny,ii)).*(mtx_P_y_tilde(:,(ii-1)*Nz+1:ii*Nz).');
        end

        vec_z((i-1)*N+1:i*N) = reshape(fft(vec_z_temp),Nx*Ny*Nz,1) / sqrt(N);
    end
end