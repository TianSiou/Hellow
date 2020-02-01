function B = FAME_Matrix_B_Anisotropic_fem_2( B, C_1, C_2, C_3, delta_1, delta_2, delta_3 )
    N = length(C_1);

    I = speye(N);

    B.B_11 = 0.5*( B.B_eps_11 + B.B_eps_11_shift_x );
    B.B_22 = 0.5*( B.B_eps_22 + B.B_eps_22_shift_y );
    B.B_33 = 0.5*( B.B_eps_33 + B.B_eps_33_shift_z );
    
    B.B_12 = 0.25*( B.B_eps_12*(2*I + delta_2*C_2') + B.B_eps_12_shift_x*(2*I + delta_2*C_2')*(I+delta_1*C_1) );
    B.B_13 = 0.25*( B.B_eps_13*(2*I + delta_3*C_3') + B.B_eps_13_shift_x*(2*I + delta_3*C_3')*(I+delta_1*C_1) );
    
    B.B_21 = 0.25*( B.B_eps_21*(2*I + delta_1*C_1') + B.B_eps_21_shift_y*(2*I + delta_1*C_1')*(I+delta_2*C_2) );
    B.B_23 = 0.25*( B.B_eps_23*(2*I + delta_3*C_3') + B.B_eps_23_shift_y*(2*I + delta_3*C_3')*(I+delta_2*C_2) );
    
    B.B_31 = 0.25*( B.B_eps_31*(2*I + delta_1*C_1') + B.B_eps_31_shift_z*(2*I + delta_1*C_1')*(I+delta_3*C_3) );
    B.B_32 = 0.25*( B.B_eps_32*(2*I + delta_2*C_2') + B.B_eps_32_shift_z*(2*I + delta_2*C_2')*(I+delta_3*C_3) );
    
%     B.lssvr = 'chol_amd';
    B.lssvr = 'pcg';
    switch B.lssvr
        case 'lu'
            B.B    = [ real(B.B_11) B.B_12       B.B_13       ;
                       B.B_21       real(B.B_22) B.B_23       ;
                       B.B_31       B.B_32       real(B.B_33)];
            [ B.L, B.U ] = lu(B.B);  
        case 'chol'
            B.B    = [ real(B.B_11) B.B_12       B.B_13       ;
                       B.B_21       real(B.B_22) B.B_23       ;
                       B.B_31       B.B_32       real(B.B_33)];
            B.U = chol(B.B);
            B.L = B.U';
        case 'lu_amd'
            B.B    = [ real(B.B_11) B.B_12       B.B_13       ;
                       B.B_21       real(B.B_22) B.B_23       ;
                       B.B_31       B.B_32       real(B.B_33)];
            B.P             = amd(B.B);
            B.invP(B.P) = 1:3*N;
            [ B.L, B.U ] = lu(B.B(B.P,B.P));  
        case 'chol_amd'
            B.B    = [ real(B.B_11) B.B_12       B.B_13       ;
                       B.B_21       real(B.B_22) B.B_23       ;
                       B.B_31       B.B_32       real(B.B_33)];
            B.P             = amd(B.B);
            B.invP(B.P) = 1:3*N;
            B.L             = chol(B.B(B.P,B.P),'lower');
            B.U             = B.L';
        case 'lumping_1'
            B.B    = [ real(B.B_11) B.B_12       B.B_13       ;
                       B.B_21       real(B.B_22) B.B_23       ;
                       B.B_31       B.B_32       real(B.B_33)];
            B.B = diag(sum(B.B.'));
            [ B.L, B.U ] = lu(B.B);  
        case 'lumping_2'
            B.B_12 = diag( sum(B.B_12.') );
            B.B_13 = diag( sum(B.B_13.') );
            B.B_21 = diag( sum(B.B_21) );
            B.B_23 = diag( sum(B.B_23.') );
            B.B_31 = diag( sum(B.B_31) );
            B.B_32 = diag( sum(B.B_32) );

            B.B    = [ real(B.B_11) B.B_12       B.B_13       ;
                       B.B_21       real(B.B_22) B.B_23       ;
                       B.B_31       B.B_32       real(B.B_33)];

            [ B.L, B.U ] = lu(B.B); 
        otherwise
%             B.B    = [ real(B.B_11) B.B_12       B.B_13       ;
%                        B.B_21       real(B.B_22) B.B_23       ;
%                        B.B_31       B.B_32       real(B.B_33)];
            B.B    = [ B.B_11       B.B_12       B.B_13       ;
                       B.B_21       B.B_22       B.B_23       ;
                       B.B_31       B.B_32       B.B_33      ];
    end
    
    if norm(B.B-B.B','fro') < 1e-12
        B.flag_hermitian = 1;
    else
        B.flag_hermitian = 0;
    end
    
end