function vec_y = FAME_Matrix_Vector_Production_FFT_Single_General(vec_x,Nx,Ny,Nz,N,D_kx,D_ky,D_kz)
    for ii = 1 : 1
        vec_x_tmp = vec_x(   (ii-1)*N + 1 : ii*N, 1 );
        vec_y_tmp = zeros(length(vec_x_tmp),1);    
        mtx_G_hat = zeros(Ny,Nz,Nx);
        for mm= 1:Nz
            mtx_vec_x         = reshape(vec_x_tmp((mm-1)*Nx*Ny+1:mm*Nx*Ny), Nx, Ny);
            mtx_tmp           = sparse(1:Nx,1:Nx,conj(D_kx)) * mtx_vec_x;
            mtx_FFT           = fft(mtx_tmp);
            mtx_G_hat(:,mm,:) = mtx_FFT';
        end
        for jj = 1:Nx
            mtx_tmp     = conj(sparse(1:Ny,1:Ny,D_ky(:,jj)) * mtx_G_hat(:,:,jj));
            mtx_FFT     = fft(mtx_tmp);
            mtx_G_tilde = mtx_FFT';
            mtx_tmp     = conj(D_kz(:,1:Ny,jj).*mtx_G_tilde);
            mtx_G       = fft(mtx_tmp);                     

            vec_y_tmp((jj-1)*Ny*Nz+1:jj*Ny*Nz,1) = reshape(mtx_G,Ny*Nz,1);
        end
        vec_y(   (ii-1)*N + 1 : ii*N, 1 ) = vec_y_tmp;
    end
    vec_y = vec_y/sqrt(N) ;
end