function vec_y = FAME_Matrix_Vector_Production_IFFT_Single_Simple(vec_x,Nx,Ny,Nz,N,D_k)

        vec_x = vec_x(       1 :   N, 1 );
        % Reshape the input vector to do the inverse fast Fourier transform
        mat_x = reshape(vec_x,Nx,Ny,Nz);
        % Inverse Fast Fourier Transform 
        vec_y = reshape(ifftn(mat_x),N,1);
        %
        vec_y = D_k.*vec_y;
        % Assemble output vector
        vec_y = sqrt(N)*vec_y;            
end