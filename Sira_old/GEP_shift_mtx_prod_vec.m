function y = GEP_shift_mtx_prod_vec(vec, A, B, sigma )

if isa(A,'function_handle')
    MtxVecMult_A = A;
else
    % Check matrix and right hand side vector inputs have appropriate sizes
    [m,n] = size(A);
    if (m ~= n)
        fprintf('error:NonSquareMatrix\n');return;
    end
    MtxVecMult_A = @(x)Amtrix(x,A); 
end

if isa(B,'function_handle')
    MtxVecMult_B = B;
else
    % Check matrix and right hand side vector inputs have appropriate sizes
    [m,n] = size(B);
    if (m ~= n)
        fprintf('error:NonSquareMatrix\n');return;
    end
    MtxVecMult_B = @(x)Bmtrix(x,B); 
end

y = MtxVecMult_A(vec) - sigma * MtxVecMult_B(vec);

end

% =======================================
%
 function u = Amtrix(x, A)
   u = A * x;
 end
 
 function u = Bmtrix(x, B)
   u = B * x;
 end
