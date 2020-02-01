function Err = FAME_Error_Check_Anisotropic( Freq_array, Ele_field_mtx,C, Cs, B )    
    Err = zeros(size(Ele_field_mtx,2),1);
    N   = size(C,1)/3;
    for i = 1:size(Ele_field_mtx,2)        
        Ele_field = Ele_field_mtx(:,i);
        
%         temp = [ B_eps.d11.*Ele_field(1:N) + B_eps.d12.*Ele_field(N+1:2*N) + B_eps.d13.*Ele_field(2*N+1:3*N);
%                  B_eps.d21.*Ele_field(1:N) + B_eps.d22.*Ele_field(N+1:2*N) + B_eps.d23.*Ele_field(2*N+1:3*N);
%                  B_eps.d31.*Ele_field(1:N) + B_eps.d32.*Ele_field(N+1:2*N) + B_eps.d33.*Ele_field(2*N+1:3*N) ];
        temp = FAME_Matrix_Vector_Production_B_Anisotropic(Ele_field, B, N);

        Err(i) = norm( Cs*(C*Ele_field) - Freq_array(i)^2*temp);
    end
end