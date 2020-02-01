function Err = FAME_Error_Check_Biisotropic_fem( Freq_array, Ele_field_mtx, Mag_field_mtx, C, Cs, C_x, C_y, C_z, B )    
    Err = zeros(size(Ele_field_mtx,2),1);
    C_w = blkdiag(C_x, C_y, C_z);
    flag_trapezoidal_approx = 0;
    for i = 1:size(Ele_field_mtx,2)
        Ele_field = Ele_field_mtx(:,i);
        Mag_field = Mag_field_mtx(:,i);
        if flag_trapezoidal_approx == 0
            Err(i) = norm( [C_w*C*Ele_field; Cs*C_w*Mag_field] - Freq_array(i)*1i*[B.Zeta*Ele_field + B.Mu*Mag_field; -B.Eps*Ele_field - B.Xi*Mag_field ] );
        else
            Err(i) = norm( [C*Ele_field; Cs*Mag_field] - Freq_array(i)*1i*[B.Zeta*Ele_field + B.Mu*Mag_field; -B.Eps*Ele_field - B.Xi*Mag_field ] );
        end
            
    end
end