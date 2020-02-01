function vec_y = FAME_Matrix_Vector_Production_invB_Anisotropic(vec_x, B, N)
    vec_y = [ B.invB_eps.d11.*vec_x(1:N) + B.invB_eps.d12.*vec_x(N+1:2*N) + B.invB_eps.d13.*vec_x(2*N+1:3*N);
              B.invB_eps.d21.*vec_x(1:N) + B.invB_eps.d22.*vec_x(N+1:2*N) + B.invB_eps.d23.*vec_x(2*N+1:3*N);
              B.invB_eps.d31.*vec_x(1:N) + B.invB_eps.d32.*vec_x(N+1:2*N) + B.invB_eps.d33.*vec_x(2*N+1:3*N) ];
end