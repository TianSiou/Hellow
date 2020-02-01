function vec_y = FAME_Matrix_Vector_Production_FFT_Triple_Simple(vec_x,Nx,Ny,Nz,N,D_ks)

        vec_x_1 = D_ks.*vec_x(         1 :   N, 1 );
        vec_x_2 = D_ks.*vec_x(     N + 1 : 2*N, 1 );
        vec_x_3 = D_ks.*vec_x(   2*N + 1 : 3*N, 1 );
        % Reshape the input vector to do the fast Fourier transform
        mat_x_1 = reshape(vec_x_1,Nx,Ny,Nz);
        mat_x_2 = reshape(vec_x_2,Nx,Ny,Nz);
        mat_x_3 = reshape(vec_x_3,Nx,Ny,Nz);
        % Fast Fourier Transform
        vec_y_1 = reshape(fftn(mat_x_1),N,1);
        vec_y_2 = reshape(fftn(mat_x_2),N,1);
        vec_y_3 = reshape(fftn(mat_x_3),N,1);
        % Assemble output vector
        vec_y = (1/sqrt(N))*[vec_y_1;vec_y_2;vec_y_3];   

end