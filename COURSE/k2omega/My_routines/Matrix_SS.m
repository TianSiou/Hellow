function mtx_ss = Matrix_SS( mtx_symplecitc, N )
    O_3N = sparse(3*N,3*N);
    I_3N = speye(3*N);
    J    = [O_3N,I_3N;-I_3N,O_3N];

    mtx_ss.fun_Ms = @(freq) mtx_symplecitc.fun_M(freq)*J*mtx_symplecitc.fun_L(freq).' + mtx_symplecitc.fun_L(freq)*J*mtx_symplecitc.fun_M(freq).';
    mtx_ss.fun_Ls = @(freq) mtx_symplecitc.fun_L(freq)*J*mtx_symplecitc.fun_L(freq).';
end