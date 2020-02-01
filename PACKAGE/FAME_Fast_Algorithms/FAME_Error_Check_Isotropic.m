function Err = FAME_Error_Check_Isotropic( Freq_array, Ele_field_mtx,C, Cs, B )    
    Err = zeros(size(Ele_field_mtx,2),1);
    for i = 1:size(Ele_field_mtx,2)        
        Ele_field = Ele_field_mtx(:,i);
%         Err(i) = norm( Cs*(C*Ele_field) - Freq_array(i)^2*B_eps.*Ele_field );
        Err(i) = norm( Cs*(C*Ele_field) - Freq_array(i)^2*FAME_Matrix_Vector_Production_B_Isotropic(Ele_field, B) );
        
    end
end