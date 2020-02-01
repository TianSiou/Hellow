function mtx_SHH = Matrix_SHH(mtx_gyroscopic)
    N = length(mtx_gyroscopic.G);

    O_3N = sparse(N,N);
    I_3N = speye(N);
    
    mtx_SHH.fun_H    = @(freq) [                      O_3N,   -mtx_gyroscopic.K; 
                                mtx_gyroscopic.fun_M(freq),                O_3N];
    mtx_SHH.fun_S    = @(freq) [mtx_gyroscopic.fun_M(freq),           mtx_gyroscopic.G;
                                                      O_3N, mtx_gyroscopic.fun_M(freq)];
    mtx_SHH.fun_S_1  = @(freq) [ I_3N,        .5*mtx_gyroscopic.G;
                                 O_3N, mtx_gyroscopic.fun_M(freq)];
    mtx_SHH.fun_S_2  = @(freq) [mtx_gyroscopic.fun_M(freq), .5*mtx_gyroscopic.G;
                                                      O_3N,                I_3N];
    
end