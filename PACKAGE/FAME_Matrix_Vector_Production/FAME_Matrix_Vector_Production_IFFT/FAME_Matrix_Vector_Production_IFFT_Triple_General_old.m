function vec_y = FAME_Matrix_Vector_Production_IFFT_Triple_General(vec_x,Nx,Ny,Nz,N,D_kx,D_ky,D_kz)
    for ii = 1 : 3    
        vec_x_tmp     = vec_x(   (ii-1)*N + 1 : ii*N, 1 );
        vec_y_tmp     = zeros(length(vec_x_tmp),1);    
        mtx_Q_hat     = zeros(Ny,Nx,Nz);
        for jj = 1:Nx
            mtx_vec_x   = reshape(vec_x_tmp((jj-1)*Ny*Nz+1:jj*Ny*Nz), Nz, Ny);
            mtx_Q_tilde = ifft(mtx_vec_x);
            mtx_tmp     = mtx_Q_tilde;

            mtx_Q_tilde       = D_kz(:,:,jj).*mtx_tmp;
            mtx_Q_tilde       = mtx_Q_tilde.';

            mtx_tmp           = ifft(mtx_Q_tilde);

            mtx_Q_hat(:,jj,:) = sparse(1:Ny,1:Ny,D_ky(:,jj)) * mtx_tmp; 
        end   
        for mm= 1:Nz
            mtx_tmp = mtx_Q_hat(:,:,mm);
            mtx_FFT = ifft(mtx_tmp.');
            mtx_Q   = sparse(1:Nx,1:Nx,D_kx) * mtx_FFT;  

            vec_y_tmp((mm-1)*Nx*Ny+1:mm*Nx*Ny,1) = reshape(mtx_Q,Nx*Ny,1);
        end
        vec_y(   (ii-1)*N + 1 : ii*N, 1 ) = vec_y_tmp;
    end
    vec_y = vec_y*sqrt(N);
end