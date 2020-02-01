function vec_y = FAME_Matrix_Vector_Production_FFT_Single_Simple(vec_x,Nx,Ny,Nz,N,D_ks)

        vec_x = D_ks.*vec_x(         1 :   N, 1 );
        % Reshape the input vector to do the fast Fourier transform
        mat_x = reshape(vec_x,Nx,Ny,Nz);
        % Fast Fourier Transform
        vec_y = reshape(fftn(mat_x),N,1);
        % Assemble output vector
        vec_y = (1/sqrt(N))*vec_y;   
end