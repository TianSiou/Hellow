
grid_num       = [8,8,8];
wave_vec       = [0,0,0.25];
k              = 2; %%number of ncell
shift          = 3;
Nx = grid_num(1); Ny = grid_num(2); Nz = grid_num(3);
N = Nx*Ny*Nz;
ps         = permu_uni2S(k, 3*k*N); %%%x(ps) is px
p          = permu_S2uni(k, 3*k*N); %%%x(p)  is p'x

B_1     = mtx_16.B.B_eps;     B_2    = mtx_9.B.B_eps;
B_eps_p   = [];
Sigma_r_u = mtx_16.Lambdas.Sigma_r;
Sigma_r   = []; 
for i = 1 : k
    if i <= k/2
        B = B_1;
    else
        B = B_2;
    end
    B_eps_p = [B_eps_p; B];
    Sigma_r = [Sigma_r; Sigma_r_u];
end
    B_eps    = B_eps_p(ps);
    fun_Ar        = @(x) FAME_Matrix_Vector_Production_Isotropic_Ar_Simple(x, B_eps, Nx, Ny, k*Nz, k*N, mtx_sup.Lambdas.Pi_Qr, mtx_sup.Lambdas.Pi_Qrs, mtx_sup.Lambdas.D_k, mtx_sup.Lambdas.D_ks);
    fun_A         = @(x) FAME_Matrix_Vector_Production_Isotropic_shift_A_Simple(x, fun_Ar, mtx_sup.Lambdas.Sigma_r, shift);
    
    fun_Ar_u1     = @(x) FAME_Matrix_Vector_Production_Isotropic_Ar_Simple(x, B_1, Nx, Ny, Nz, N, mtx_16.Lambdas.Pi_Qr, mtx_16.Lambdas.Pi_Qrs, mtx_16.Lambdas.D_k, mtx_16.Lambdas.D_ks);
    fun_A_u1      = @(x) FAME_Matrix_Vector_Production_Isotropic_shift_A_Simple(x, fun_Ar_u1, mtx_16.Lambdas.Sigma_r, shift);
    fun_Ar_u2     = @(x) FAME_Matrix_Vector_Production_Isotropic_Ar_Simple(x, B_2, Nx, Ny, Nz, N, mtx_9.Lambdas.Pi_Qr, mtx_9.Lambdas.Pi_Qrs, mtx_9.Lambdas.D_k, mtx_9.Lambdas.D_ks);
    fun_A_u2      = @(x) FAME_Matrix_Vector_Production_Isotropic_shift_A_Simple(x, fun_Ar_u2, mtx_16.Lambdas.Sigma_r, shift);
    fun_invA_u1   = @(x) FAME_Matrix_Vector_Production_Isotropic_shift_invA_Simple(x, fun_Ar_u1, mtx_16.Lambdas.Sigma_r, shift); 
    fun_invA_u2   = @(x) FAME_Matrix_Vector_Production_Isotropic_shift_invA_Simple(x, fun_Ar_u2, mtx_9.Lambdas.Sigma_r, shift);

    fun_P_tile   = @(x) permu_uni2S_tile(x, Nx, Ny, Nz, N, mtx_16.Lambdas, B_eps, B_1, B_2, Sigma_r, k, ps, p);
    fun_Ps_tile  = @(x) permu_S2uni_tile(x, Nx, Ny, Nz, N, mtx_16.Lambdas, B_eps, B_1, B_2, Sigma_r, k, ps, p);
    fun_M        = @(x) Matrix_Vector_Production_M_Simple(x, Nx, Ny, Nz, N, mtx_16.Lambdas, B_eps, Sigma_r, k, p, ps) - shift *x;
    fun_inv_M    = @(x) Matrix_Vector_Production_invM_Simple(x, N, fun_invA_u1, fun_invA_u2, fun_P_tile,fun_Ps_tile, k);
%     fun_M_test   = @(x) Matrix_Vector_Production_invM_Simple_test(x, N, fun_A_u1, fun_A_u2, fun_P_tile, fun_Ps_tile, k);
    fun_inv_M_test1  = @(x)  Matrix_Vector_Production_invM_Simple_test1(x, N, fun_Ar_u1, fun_Ar_u2, mtx_16.Lambdas.Sigma_r, k, shift);
%     fun_A_tile   = @(x) (fun_inv_M( fun_A(x) ) );%%%normal
    fun_A_tile   = @(x) fun_inv_M_test1( fun_Ar(x) - shift * ( (mtx_sup.Lambdas.Sigma_r.^2) .\ x ) );
    b            = rand(2*N*k,1);        
%     b_tile       = fun_inv_M(b);%%%normal
%     [ y, flag, res, iter, resvec] = bicgstabl(fun_A_tile, b_tile, 1e-12, 200); %%%normal
    y       = mtx_sup.Lambdas.Sigma_r .\ b;
    [ yy, flag, res, iter, resvec] = bicgstabl(fun_A_tile, fun_inv_M_test1(y), 1e-12, 1000);
    y       = mtx_sup.Lambdas.Sigma_r .\ yy;
%     [ y, flag, res, iter, resvec] = minres(fun_A, b, 1e-12, 200, fun_inv_M);
%     [ y, flag, res, iter, resvec] = minres(fun_A, b, 1e-12, 200);
    r = norm( fun_A(y) - b);
 function vec_x = Matrix_Vector_Production_invM_Simple_test1(vec_x, N, fun_A_u1, fun_A_u2, Sigma_r, k, shift)
    for i = 1 : k
        y = vec_x(1+(i-1)*2*N : i*2*N);
     if i <= k/2
         vec_x(1+(i-1)*2*N : i*2*N) = minres(@(x) fun_A_u1(x) - shift* ((Sigma_r.^2) .\ x), y, 1e-10, 1000);
     else
         vec_x(1+(i-1)*2*N : i*2*N) = minres(@(x) fun_A_u2(x) - shift* ((Sigma_r.^2) .\ x), y, 1e-10, 1000);
     end
    end
end  
function vec_y = Matrix_Vector_Production_invM_Simple(x, N, fun_invA_u1, fun_invA_u2, fun_P_tile, fun_Ps_tile, k)
    vec_y = fun_Ps_tile(x);
    for i = 1 : k
     if i <= k/2
         vec_y(1+(i-1)*2*N : i*2*N) = fun_invA_u1( vec_y(1+(i-1)*2*N : i*2*N) );
     else
         vec_y(1+(i-1)*2*N : i*2*N) = fun_invA_u2( vec_y(1+(i-1)*2*N : i*2*N) );
     end
    end
    vec_y = fun_P_tile(vec_y);
end
function vec_y = FAME_Matrix_Vector_Production_Isotropic_shift_invA_Simple(vec_x, fun_Ar, Sigma_r, shift) 
    vec_y = Sigma_r .\ vec_x;
    vec_y = minres(@(x) fun_Ar(x) - shift * ( (Sigma_r.^2) .\ x), vec_y, 1e-12, 1000);    
    vec_y = Sigma_r .\ vec_y;
end
function vec_y = FAME_Matrix_Vector_Production_Isotropic_shift_A_Simple(vec_x, fun_Ar, Sigma_r, shift) %% sima*QBQ*sigma
    vec_y = Sigma_r .* vec_x;
    vec_y = fun_Ar(vec_y) - shift * ( (Sigma_r.^2) .\ vec_y);
    vec_y = Sigma_r .* vec_y;
end
function vec_y = FAME_Matrix_Vector_Production_Isotropic_Ar_Simple(vec_x, B, Nx, Ny, Nz, N, Pi_Qr, Pi_Qrs, D_k, D_ks) %%%QBQ
%     vec_y = Sigma_r .* vec_x;
    vec_y = FAME_Matrix_Vector_Production_Qr_Simple(vec_x, Nx, Ny, Nz, N, Pi_Qr, Pi_Qrs, D_k, D_ks, 'normal');
  
    vec_y = B .\ vec_y;

    vec_y = FAME_Matrix_Vector_Production_Qr_Simple(vec_y, Nx, Ny, Nz, N, Pi_Qr, Pi_Qrs, D_k, D_ks, 'hermitian');
%     vec_y = Sigma_r .* vec_y;
end
function vec_z = Matrix_Vector_Production_M_Simple(x, Nx, Ny, Nz, N, Lambdas, B_eps, Sigma_r, k, p, ps)
N2    = 2*N;
N3    = 3*N;
vec_y = [];
vec_z = [];
x = Sigma_r .* x;
 for i = 1 : k
    temp = FAME_Matrix_Vector_Production_Qr_Simple(x(1+(i-1)*N2 : i*N2), Nx, Ny, Nz, N, Lambdas.Pi_Qr, Lambdas.Pi_Qrs, Lambdas.D_k, Lambdas.D_ks, 'normal');
    vec_y = [vec_y; temp];
 end
    vec_y = vec_y(ps);
    vec_y = B_eps .\ vec_y;
    vec_y = vec_y(p);
for i = 1 : k
    temp = FAME_Matrix_Vector_Production_Qr_Simple(vec_y(1+(i-1)*N3 : i*N3), Nx, Ny, Nz, N, Lambdas.Pi_Qr, Lambdas.Pi_Qrs, Lambdas.D_k, Lambdas.D_ks, 'hermitian');
    vec_z = [vec_z; temp];
end
 vec_z = Sigma_r .* vec_z;
end

function vec_z = permu_uni2S_tile(x, Nx, Ny, Nz, N, Lambdas, B_eps, B_1, B_2, Sigma_r, k, ps, p) %%%P_tile
N2  = 2*N;
N3  = 3*N;
vec_y = [];
vec_z = [];
for i = 1 : k
 if i <= k/2
     temp = x(1+(i-1)*N2 : i*N2) .* Lambdas.Sigma_r;
     temp = FAME_Matrix_Vector_Production_Qr_Simple(temp, Nx, Ny, Nz, N, Lambdas.Pi_Qr, Lambdas.Pi_Qrs, Lambdas.D_k, Lambdas.D_ks, 'normal');
     temp = B_1 .\ temp;
     vec_y = [vec_y; temp];
 else
     temp = x(1+(i-1)*N2 : i*N2) .* Lambdas.Sigma_r;
     temp = FAME_Matrix_Vector_Production_Qr_Simple(temp, Nx, Ny, Nz, N, Lambdas.Pi_Qr, Lambdas.Pi_Qrs, Lambdas.D_k, Lambdas.D_ks, 'normal');
     temp = B_2 .\ temp;   
     vec_y = [vec_y; temp];
 end
end
    vec_y = vec_y(ps);
    vec_y = B_eps .* vec_y;
    vec_y = vec_y(p);
    for i = 1 : k
        temp = FAME_Matrix_Vector_Production_Qr_Simple(vec_y(1+(i-1)*N3 : i*N3), Nx, Ny, Nz, N, Lambdas.Pi_Qr, Lambdas.Pi_Qrs, Lambdas.D_k, Lambdas.D_ks, 'hermitian');
        vec_z = [vec_z; temp];
    end
    vec_z = Sigma_r .\ vec_z;
end

function vec_z = permu_S2uni_tile(x, Nx, Ny, Nz, N, Lambdas, B_eps, B_1, B_2, Sigma_r, k, ps, p) %%%P'_tile
N2  = 2*N;
N3  = 3*N;
vec_y = [];
vec_z = [];
y = Sigma_r .\ x;
for i = 1 : k
    temp  = FAME_Matrix_Vector_Production_Qr_Simple( y(1+(i-1)*N2 : i*N2), Nx, Ny, Nz, N, Lambdas.Pi_Qr, Lambdas.Pi_Qrs, Lambdas.D_k, Lambdas.D_ks, 'normal');
    vec_y = [vec_y; temp]; 
end
vec_y = vec_y(ps);
vec_y = B_eps .* vec_y;
vec_y = vec_y(p);
for i = 1 : k
    if i <= k/2
         temp  = B_1 .\ vec_y(1 + (i-1)*N3 : i*N3);
         temp  = FAME_Matrix_Vector_Production_Qr_Simple(temp, Nx, Ny, Nz, N, Lambdas.Pi_Qr, Lambdas.Pi_Qrs, Lambdas.D_k, Lambdas.D_ks, 'hermitian');
         temp  = Lambdas.Sigma_r .* temp;
         vec_z = [vec_z; temp];
    else
         temp = B_2 .\ vec_y(1 + (i-1)*N3 : i*N3);
         temp = FAME_Matrix_Vector_Production_Qr_Simple(temp, Nx, Ny, Nz, N, Lambdas.Pi_Qr, Lambdas.Pi_Qrs, Lambdas.D_k, Lambdas.D_ks, 'hermitian');
         temp = Lambdas.Sigma_r .* temp;      
         vec_z = [vec_z; temp];
    end
end
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
    
    
function vec_y = permu_S2uni(ncell, N)    %%%P'
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