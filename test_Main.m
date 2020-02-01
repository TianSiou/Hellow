
    idouble_pi = 1i*2*pi;
    grid_num       = [8,8,8];
    wave_vec       = [0,0,0.25];
    k              = 2; %%number of ncell
    lattice_vec_a3 = [0, 0 ,1];
    shift          = 0.1;
    Nx = grid_num(1); Ny = grid_num(2); Nz = grid_num(3);
    one_z      = ones(k * grid_num(3),1);
    
    K_z        = spdiags([-one_z one_z],0:1,k * grid_num(3),k * grid_num(3));
    theta_z    = dot(lattice_vec_a3,wave_vec);
    K_z(end,1) = exp(idouble_pi * theta_z * k);
    hz = 1/8;
    K_z = K_z/hz;
    
    N = Nx*Ny*Nz;
    if norm(wave_vec) == 0
        Nd = 2*N - 2;
    else 
        Nd = 2*N;
    end
    O = sparse(k * N, k * N);
    C_1     = mtx_16.C_1; 
    C_2     = mtx_16.C_2;  
    C_3     = mtx_16.C_3;
    C_u     = mtx_16.C; 
    Cs_u    = C_u';
    C = [O,                        -( kron(K_z, kron(speye(Ny), speye(Nx) ) ) ),     kron(speye(k), C_2);
         kron(K_z, kron(speye(Ny), speye(Nx) ) ),         O,                     -(kron(speye(k), C_1) ) ;
         -( kron(speye(k), C_2) ),                  kron(speye(k), C_1),                       O];
% C = [O,                        -( kron(kron(speye(Ny), speye(Nx) ), K_z ) ),     kron(speye(k), C_2);
%          kron(kron(speye(Ny), speye(Nx) ), K_z ),         O,                     -(kron(speye(k), C_1) ) ;
%          -( kron(speye(k), C_2) ),                  kron(speye(k), C_1),                       O];
    Cs = C';  
    C_M     = [O , -(kron(speye(k), C_3) ) , kron(speye(k) , C_2);
               kron(speye(k), C_3) , O , -(kron(speye(k), C_1) );
               -(kron(speye(k), C_2) ) , kron(speye(k) , C_1), O];
    Cs_M    = C_M';
    
   p         = ( kron(ones(3*k,1), [1:N]' )) + N *kron(  kron([1:k]' -1,ones(3,1) ) + k * kron( ones(k,1), [1:3]' -1) ,ones(N,1)) ;  
   ps         = permu_uni2S(k, 3*N*k); %%% uni to sup
    %%%%test CM
%     x       = randn(3*ncell*N);
    C_p    = kron(speye(k) , C_u);
%     C_p     = Ps(test, ncell, 3*ncell*N, 3*ncell*N);
%     C_p     = P(C_p, ncell, 3*ncell*N, 3*ncell*N);
%       z       = norm(full(C_p(ps,ps))-full(C_M));
    B_1     = mtx_16.B.B_eps;     B_2    = mtx_9.B.B_eps;
    B_eps_p   = [];
    for i = 1 : k
        if i <= k/2
            B = B_1;
        else
            B = B_2;
        end
        B_eps_p = [B_eps_p; B];
    end
%     B_esp    = P(B_esp_p, ncell, 3*ncell*N, 1);
    B_eps    = B_eps_p(ps);
    B_eps    = spdiags(B_eps, 0 , 3*k*N,3*k*N);
    alpha_1  = mean(B_1);
    alpha_2  = mean(B_2);
    alpha_Beps = mean(B_eps);
    B_p      = spdiags(B_eps_p, 0 , k*3*N,k*3*N);
    B_1      = spdiags(B_1, 0 , 3*N, 3*N);
    B_2      = spdiags(B_2, 0 , 3*N, 3*N);   
    M        = Cs_M*C_M - shift * B_eps;
%       M        = Cs*C - alpha_Beps * speye(3*k*N);
    %%%test M
%    x         = rand(3*ncell*N,1);
%    test = C_p'*C_p-shift*B_p;
%   test = ( kron(speye(ncell) , C_u)'*kron(speye(ncell) , C_u) - shift*B_p);
%   y    = x(p);
%   y    = test * y;
%   y    = y(ps);
%   z    = norm(M *x - y);
%   z    = norm(M *x - test(ps,ps)*x);

    Sigma_1  = spdiags(mtx_16.Lambdas.Sigma, 0 , 3*N, 3*N);
    Sigma_2  = spdiags(mtx_9.Lambdas.Sigma, 0 , 3*N, 3*N);

    

%    fun_A           =  @(x) (Cs * C - shift * B_eps) * x;
   fun_A           =  @(x) (mtx_sup.Cs * mtx_sup.C - shift * B_eps )*x;
   fun_A_u1        =  @(x) (Cs_u * C_u - shift * B_1) *x;
   fun_A_u2        =  @(x) (Cs_u * C_u - shift * B_2) *x;
   fun_preA_u1     =  @(x) Matrix_Vector_Production_preCu_Simple(x, Nx, Ny, Nz, N, mtx_16.Lambdas, alpha_1);
   fun_preA_u2     =  @(x) Matrix_Vector_Production_preCu_Simple(x, Nx, Ny, Nz, N, mtx_9.Lambdas , alpha_2);
   fun_A_u1_tile   =  @(x) fun_preA_u1( fun_A_u1(x));
   fun_A_u2_tile   =  @(x) fun_preA_u2( fun_A_u2(x));
   

   
   fun_invA_u1     =  @(x) bicgstabl(fun_A_u1_tile , fun_preA_u1(x), 1e-12, 1000);
   fun_invA_u2     =  @(x) bicgstabl(fun_A_u2_tile , fun_preA_u2(x), 1e-12, 1000);
   %%%%test fun_invA_u1&u2
%    x = rand(3*N,1);
%    xx = rand(3*N,1);
%    y = fun_invA_u1(x);
%    yy = fun_invA_u2(xx); 
%    z = norm(fun_A_u1_tile(y) - fun_preA_u1(x) );
%    zz = norm(fun_A_u2_tile(yy) - fun_preA_u2(xx) );
   
   
   fun_invM        =  @(x) Matrix_Vector_Production_invM_Simple(x, fun_invA_u1, fun_invA_u2, k, N, ps, p);
%    fun_invM        = @(x)  Matrix_Vector_Production_preC_Simple(x, Nx, Ny, k*Nz, k*N, mtx_sup.Lambdas.Pi_Q, mtx_sup.Lambdas.Pi_Qs, mtx_sup.Lambdas.D_k, mtx_sup.Lambdas.D_ks, mtx_sup.Lambdas.Sigma, alpha_Beps);
   %%%test invM
   x = rand(3*N*k,1);
   x = x / norm(x);
   y = fun_invM(x);
   z = norm(M*y - x);%%再造一個super M
   fun_A_tile      =  @(x) fun_invM( fun_A(x) );
   b               = rand(3*k*N,1);
   b_tile          =  fun_invM(b);
   tic
%    [ y, flag, res, iter, resvec] = minres(fun_A, b, 1e-12, 1000);
%    [ y, flag, res, iter, resvec] = minres(fun_A, b,1e-12, 100, fun_invM);
   
   [ y, flag, res, iter, resvec] = bicgstabl(fun_A_tile, b_tile, 1e-8, 100);
%    [ y, flag, res, iter, resvec] = minres(fun_A, b,1e-12, 100, fun_invM);
   toc
    r              = norm( fun_A(y) - b )
 
%    [ y, flag, res, iter, resvec] = bicgstabl(@(x)speye(3*k*N)* x - Matrix_Vector_Production_QsBQ_Simple(x, Nx, Ny, k*Nz, k*N, B_eps_p(ps), mtx_sup.Lambdas, shift), b, 1e-8, 100);
    
    
    
% function vec_y = Matrix_Vector_Production_QsBQ_Simple(vec_x, Nx, Ny, Nz, N, B, Lambdas, shift)
%     vec_y  = Lambdas.Sigma .\ vec_x;
%     vec_y  =  FAME_Matrix_Vector_Production_Q_Simple(vec_y,Nx,Ny,Nz,N,Lambdas.Pi_Q, Lambdas.Pi_Qs, Lambdas.D_k, Lambdas.D_ks,'normal');
%     vec_y  = B .*vec_y;
%     vec_y  =  FAME_Matrix_Vector_Production_Q_Simple(vec_y,Nx,Ny,Nz,N,Lambdas.Pi_Q, Lambdas.Pi_Qs, Lambdas.D_k, Lambdas.D_ks,'hermitian');
%     vec_y  = Lambdas.Sigma .\ vec_y;
%     vec_y  = vec_y * shift;
% end
    
    
function vec_y = Matrix_Vector_Production_preC_Simple(vec_x, Nx, Ny, Nz, N, Pi_Q, Pi_Qs, D_k, D_ks, Sigma, alpha)
%     vec_y  = Qs(vec_x);
    vec_y  = FAME_Matrix_Vector_Production_Q_Simple(vec_x,Nx,Ny,Nz,N,Pi_Q,Pi_Qs,D_k,D_ks,'hermitian');
    Sigma  = spdiags(Sigma, 0 , 3*N, 3*N);
    vec_y  = (Sigma.^2 - alpha * speye(3*N) ) \ vec_y;
    vec_y  =  FAME_Matrix_Vector_Production_Q_Simple(vec_y,Nx,Ny,Nz,N,Pi_Q,Pi_Qs,D_k,D_ks,'normal');
%     vec_y  =  Q(vec_y);
end    
    
    
function vec_y = Matrix_Vector_Production_preCu_Simple(vec_x, Nx, Ny, Nz, N, Lambdas, alpha)
    vec_y  = FAME_Matrix_Vector_Production_Q_Simple(vec_x,Nx,Ny,Nz,N,Lambdas.Pi_Q,Lambdas.Pi_Qs,Lambdas.D_k,Lambdas.D_ks,'hermitian');
%     Sigma  = spdiags(Sigma, 0 , 3*N, 3*N);
    vec_y  = (Lambdas.Sigma.^2 - alpha ) .\ vec_y;

    vec_y  =  FAME_Matrix_Vector_Production_Q_Simple(vec_y,Nx,Ny,Nz,N,Lambdas.Pi_Q,Lambdas.Pi_Qs,Lambdas.D_k,Lambdas.D_ks,'normal');
end
function y = Matrix_Vector_Production_invM_Simple(vec_x, fun_invA_u1, fun_invA_u2, ncell, N, ps,p)
           y  = vec_x(p);
         partion   = (3*ncell*N) / ncell;
         for i = 1 : ncell
             if i <= ncell/2
                    y(1+ (i-1)* partion : i*partion) = fun_invA_u1( y(1+ (i-1)* partion : i*partion) );
             else
                    y(1+ (i-1)* partion : i*partion) = fun_invA_u2( y(1+ (i-1)* partion : i*partion) );
             end
         end
           y  = y(ps);
end
 
function vec_y = permu_uni2S(ncell, N)  %%%P
    idx    = zeros(N,1);
    part   = N / ncell / 3; 
    idx_B  = zeros(ncell * 3,1);
    for i = 1 : 3
        for j = 1 : ncell
            idx_B( (i-1)*ncell + j) = i  + (j-1)*3;
        end
    end
    
    for i = 1 : ncell*3
        for j = 1 : part
            idx(part*(i-1) + j) = (idx_B(i) -1)* part + j;
        end
    end
    vec_y   = idx;
end
    
    
function vec_y = permu_S2uni(x, ncell, N)    %%%P'
    idx    = zeros(N,1);
    part   = N / ncell / 3; 
    idx_B  = zeros(ncell * 3,1);
    for i = 1 : ncell
        for j = 1 : 3
            idx_B( (i-1)*3 + j) = i  + (j-1)*ncell;
        end
    end
    
    for i = 1 : ncell*3
        for j = 1 : part
            idx(part*(i-1) + j) = (idx_B(i) -1)* part + j;
        end
    end
    vec_y   = idx;
end

