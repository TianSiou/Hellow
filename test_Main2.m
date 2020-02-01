    grid_num       = [4,4,4];
    wave_vec       = [0,0,0.25];
    ncell          = 2;
    shift          = 0.5;
    Nx = grid_num(1); Ny = grid_num(2); Nz = grid_num(3);
    N = Nx*Ny*Nz;
    if norm(wave_vec) == 0
        Nd = 2*N - 2;
    else 
        Nd = 2*N;
    end
    C  = mtx_sup.C;
    Cs = C';  
   p         = ( kron(ones(3*ncell,1), [1:N]' )) + N *kron(  kron([1:ncell]' -1,ones(3,1) ) + ncell * kron( ones(ncell,1), [1:3]' -1) ,ones(N,1)) ;  
   ps         = permu_uni2S(ncell, 3*N*ncell); 
    B_1     = mtx_16.B.B_eps;     B_2    = mtx_9.B.B_eps;
    B_esp   = [];
    for i = 1 : ncell
        if i <= ncell/2
            B = B_1;
        else
            B = B_2;
        end
        B_esp = [B_esp; B];
    end
    B_esp    = B_esp(ps);
    alpha_Besp = mean(B_esp);
%     B_esp    = spdiags(B_esp, 0 , 3*ncell*N,3*ncell*N);
    M        = Cs*C - alpha_Besp * speye(3*ncell*N);

   fun_A           =  @(x) Cs * (C *x)- shift * (B_esp .* x);
   fun_invM        = @(x)  Matrix_Vector_Production_preC_Simple(x, Nx, Ny, ncell*Nz, ncell*N, mtx_sup.Lambdas.Pi_Q, mtx_sup.Lambdas.Pi_Qs, mtx_sup.Lambdas.D_k, mtx_sup.Lambdas.D_ks, mtx_sup.Lambdas.Sigma, alpha_Besp);
   fun_A_tile      =  @(x) fun_invM( fun_A(x) );

    b       = rand(2*N*ncell,1);
    MVP_Qr  = @(x) FAME_Matrix_Vector_Production_Qr_Simple(x,Nx,Ny,ncell*Nz,ncell*N,mtx_sup.Lambdas.Pi_Qr,mtx_sup.Lambdas.Pi_Qrs,mtx_sup.Lambdas.D_k,mtx_sup.Lambdas.D_ks,'normal');
    MVP_Qrs = @(x)FAME_Matrix_Vector_Production_Qr_Simple(x,Nx,Ny,ncell*Nz,ncell*N,mtx_sup.Lambdas.Pi_Qr,mtx_sup.Lambdas.Pi_Qrs,mtx_sup.Lambdas.D_k,mtx_sup.Lambdas.D_ks,'hermitian');
    b_tile  = MVP_Qr(mtx_sup.Lambdas.Sigma_r .* b);
    b_tile  =  fun_invM(b_tile);
   tic
%    [ x_tile, flag, res, iter, resvec] = bicgstabl(fun_A_tile, b_tile, 1e-12, 100);
   [ x_tile, flag, res, iter, resvec] = minres(fun_A, b,1e-12, 100, fun_invM);
   toc
   
   x        = MVP_Qrs( B_esp .* x_tile);
   x        = mtx_sup.Lambdas.Sigma_r .\ x;
   fun_Ar   = @(x)  FAME_Matrix_Vector_Production_Isotropic_Ar_Simple(x, B_esp, Nx, Ny, ncell*Nz, ncell*N, mtx_sup.Lambdas.Pi_Qr, mtx_sup.Lambdas.Pi_Qrs, mtx_sup.Lambdas.D_k, mtx_sup.Lambdas.D_ks,mtx_sup.Lambdas.Sigma_r);
    r       = norm(fun_Ar(x) - shift*x - b);
    
    
    
    
function vec_y = FAME_Matrix_Vector_Production_Isotropic_Ar_Simple(vec_x, B, Nx, Ny, Nz, N, Pi_Qr, Pi_Qrs, D_k, D_ks,Sigma)
    vec_y = Sigma .* vec_x;
    vec_y = FAME_Matrix_Vector_Production_Qr_Simple(vec_y, Nx, Ny, Nz, N, Pi_Qr, Pi_Qrs, D_k, D_ks, 'normal');

    vec_y = B .\ vec_y;

    vec_y = FAME_Matrix_Vector_Production_Qr_Simple(vec_y, Nx, Ny, Nz, N, Pi_Qr, Pi_Qrs, D_k, D_ks, 'hermitian');
    vec_y = Sigma .* vec_y;
end
   
function vec_y = Matrix_Vector_Production_preC_Simple(vec_x, Nx, Ny, Nz, N, Pi_Q, Pi_Qs, D_k, D_ks, Sigma, alpha)
%     vec_y  = Qs(vec_x);
    vec_y  = FAME_Matrix_Vector_Production_Q_Simple(vec_x,Nx,Ny,Nz,N,Pi_Q,Pi_Qs,D_k,D_ks,'hermitian');
    Sigma  = spdiags(Sigma, 0 , 3*N, 3*N);
    vec_y  = (Sigma.^2 - alpha * speye(3*N) ) \ vec_y;
    vec_y  =  FAME_Matrix_Vector_Production_Q_Simple(vec_y,Nx,Ny,Nz,N,Pi_Q,Pi_Qs,D_k,D_ks,'normal');
%     vec_y  =  Q(vec_y);
end    
    
    
function vec_y = Matrix_Vector_Production_preCu_Simple(vec_x, Nx, Ny, Nz, N, Pi_Q, Pi_Qs, D_k, D_ks, Sigma, alpha)
    vec_y  = FAME_Matrix_Vector_Production_Q_Simple(vec_x,Nx,Ny,Nz,N,Pi_Q,Pi_Qs,D_k,D_ks,'hermitian');
    Sigma  = spdiags(Sigma, 0 , 3*N, 3*N);
    vec_y  = (Sigma.^2 - alpha * speye(3*N) ) \ vec_y;

    vec_y  =  FAME_Matrix_Vector_Production_Q_Simple(vec_y,Nx,Ny,Nz,N,Pi_Q,Pi_Qs,D_k,D_ks,'normal');
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

