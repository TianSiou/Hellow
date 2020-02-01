function [z] = GT_linear_system(n, mtx_A, mtx_Q, solve_LS_PQEP, b, shift)
    %
    %  Solve linear system
    %       N z = b,
    %  where 
    %       N = ( M - shift L ) J ( M^T - shift L^T ) J^T
    %  with
    %           | 0    I |
    %       J = |        |.
    %           | -I   0 |
    %
%=    
    %
    %  (i) Solve ( M - shift L ) [ z1' z2' ]' = b
    %
    
    b1 = b(1:n) - shift*b(n+1:2*n);
    z1 = solve_LS_PQEP(b1, 'normal');
    
%     shift_value_in_GT_linear_system = shift
%     Qp = shift^2*mtx_A.' - shift*mtx_Q + mtx_A;
%     test_in_GT_linear_system = norm(Qp*z1 - b1)
%     z1 = LU_W.Perm_amd * (LU_W.upper_U \ (LU_W.Low_L \ (LU_W.Perm_LU * b1(LU_W.Perm_amd_vec,:))));
%     z2 = (Q - shift*(A.'))*z1 - b(n+1:2*n);

    if isa(mtx_Q,'function_handle')
        tmp1 = mtx_Q(z1);
    else
        tmp1 = mtx_Q * z1;
    end
    
    if isa(mtx_A,'function_handle') 
        tmp3   = mtx_A(z1, 'transp');
    else 
        tmp3 = mtx_A.' * z1;
    end  

    z2 = tmp1 - shift * tmp3 - b(n+1:2*n);
    
    %
    %  (ii) Compute [z1' z2']' := J^T [z1' z2']'
    %
    
    tmp = -z2;
    z2  = z1;
    z1  = tmp;
    
    %
    %  (iii) Solve ( M^T - shift L^T ) [ z1' z2' ]' = [ z1' z2' ]' 
    %
    
    z3 = -z2;
%     tmp = z1 - (Q - shift*A)*z3;
    
    if isa(mtx_Q,'function_handle')
        tmp1 = mtx_Q(z3);
    else
        tmp1 = mtx_Q * z3;
    end
    
    if isa(mtx_A,'function_handle')
        tmp2 = mtx_A(z3, 'notransp'); 
    else
        tmp2 = mtx_A * z3; 
    end 
    tmp = z1 - tmp1 + shift * tmp2;
    
%     z1 = LU_W.Perm_amd * (LU_W.Perm_LU_Tran * (LU_W.Low_L_Tran \ (LU_W.upper_U_Tran \ tmp(LU_W.Perm_amd_vec,:))));
    z1 = solve_LS_PQEP(tmp, 'transpose');
    z2 = z3 - shift*z1;
    
    %
    %  (iv) Compute z := J [z1' z2']'
    %
    
    z = [z2; -z1];
end