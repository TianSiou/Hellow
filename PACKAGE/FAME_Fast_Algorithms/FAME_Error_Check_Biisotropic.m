function Err = FAME_Error_Check_Biisotropic( Freq_array, Ele_field_mtx, Mag_field_mtx, C, Cs, B )    
    Err = zeros(size(Ele_field_mtx,2),1);
    for i = 1:size(Ele_field_mtx,2)
        Ele_field = Ele_field_mtx(:,i);
        Mag_field = Mag_field_mtx(:,i);
        
        Err(i) = norm( [C*Ele_field; Cs*Mag_field] - Freq_array(i)*FAME_Matrix_Vector_Production_B_Biisotropic([Ele_field;Mag_field], B) );
    end
end