
clear 
format long e;

global tol_LinearSys inner_iter_count outer_iter_count
global LU_precond flag_LU


vec_k      = zeros(3,1);
fid        = fopen('runsh_parameter.txt','r');
n1         = fscanf(fid,'%d',1); 
sigma               = fscanf(fid,'%g',1);
ncv                 = fscanf(fid,'%d',1);  
file_name           = fscanf(fid,'%s',1);
icase               = fscanf(fid,'%d',1);
icount              = fscanf(fid,'%d',1);
chirality_parameter = fscanf(fid,'%g',1);
chirality_parameter = 0;
lattice             = fscanf(fid,'%s',1);
media_type          = fscanf(fid,'%s',1);
status              = fclose(fid);

fid            = fopen( file_name, 'w+');
                      

% partitions of chiral medium parameter                 
c_size = 2; 


extmag_intensity             = 0.875;
permitt                      = 13;
permitt_2                    = permitt*extmag_intensity;     % 
permitt_1                    = sqrt(permitt^2 + permitt_2^2);
ele_permitt_in  = [      permitt_1,    0,         -1i*permitt_2;
                         0,            permitt,   0;
                         1i*permitt_2, 0,     permitt_1];
ele_permitt_out = eye(3);
mag_permeab_in  = eye(3);
mag_permeab_out = eye(3);

flag_eps = 1;

% Diagonal matries associated with dielectric function and permeability function
% chirality_parameter = 0.5; %0.8; %3.6; %chirality_parameter_array(2);

switch lattice
        
    case {'FCC', 'fcc'}
        
        cubic_length   = 1.0; %1.0 ;       %2.0d0 
        edge_len       = cubic_length / sqrt(2);

        period_vec1    = edge_len * [ 1; 0; 0 ];
        period_vec2    = edge_len * [ 1/2; sqrt(3/4); 0 ];
        period_vec3    = edge_len * [ 1/2; sqrt(1/12); sqrt(2/3) ];

        length_x       = period_vec1(1);  
        length_y       = period_vec2(2);
        length_z       = period_vec3(3);

        ratio_cyld_length   = 0.06; %0.11; %0.075; %0.11; %0.11; %0.06;
        ratio_radius_length = 0.08; %0.15; %0.12; %0.15;

        rd_sphere           = ratio_radius_length * cubic_length;
        rd_cylinder         = ratio_cyld_length   * cubic_length;
        
        if ( mod(n1,6) ~= 0 )
            n1 = (floor(n1/6) + 1) * 6;
        end
        n2             = n1; %6 * m_2;
        n3             = n1; %ceil(2 * sqrt(6) * m_1);

        delta          = zeros(3,1);     
        delta(1)       = length_x / n1;
        delta(2)       = length_y / n2;
        delta(3)       = length_z / n3;

        n              = n1 * n2 * n3;
        
        shape      = 'ellipse'; %'cylinder'; %'ellipse';
        
        [mtx_data.N, mtx_data.M] = Generate_mtx_NM_Lebedev(n1, n2, n3, edge_len, rd_sphere, rd_cylinder, ele_permitt_in, ele_permitt_out, mag_permeab_in, mag_permeab_out, delta, shape);
       
                
        no_wave_vec    = 15;
        tt             = linspace(0,1,no_wave_vec);
        
        vec_X = [0 ; 1 ; 0];
        vec_U = [ 1/4 ; 1 ; 1/4 ];
        vec_L = [ 1/2 ; 1/2 ; 1/2 ];
        vec_G = [0; 0; 0 ];
        vec_W = [ 1/2 ; 1 ; 0 ];
        vec_K = [ 3/4 ; 3/4 ; 0 ];

        mtx_Q = [ 0.5  0.5  0; -1/(2*sqrt(3))  1/(2*sqrt(3))  1/sqrt(3); 1/sqrt(6)  -1/sqrt(6)  1/sqrt(6) ];
        mtx_Q = sqrt(2) * mtx_Q;

        vec_X = mtx_Q * vec_X;
        vec_U = mtx_Q * vec_U;
        vec_L = mtx_Q * vec_L;
        vec_W = mtx_Q * vec_W;
        vec_K = mtx_Q * vec_K;
        
        switch icase
            case 1
                vec_k  =  (1-tt(icount)) * vec_X + tt(icount) * vec_U;
            case 2
                vec_k  =  (1-tt(icount)) * vec_U + tt(icount) * vec_L;
            case 3
                if icount ~= no_wave_vec 
                   vec_k  =  (1-tt(icount)) * vec_L;
                end
            case 4
                if   icount ~= 1
                    vec_k  =   tt(icount) * vec_X;
                end
            case 5
                vec_k  =  (1-tt(icount)) * vec_X + tt(icount) * vec_W;
            case 6
                vec_k  =  (1-tt(icount)) * vec_W + tt(icount) * vec_K;
        end
        
        vec_k          = vec_k / cubic_length;
        
        %
        % =============== Construct EigDecompDoubCurl_cell ========================
        % 
            
        %  Generate the matrices mtx_A and mtx_B for the generalized eigenvalue problem

        [FFT_parameter, Lambda_1, Lambda_2, Lambda_3] = construct_fft_mtx(n1, n2, n3, vec_k, period_vec1, period_vec2, period_vec3, delta);
        
        % cluster = '000'; 
        SVD_curl = SVDSingleCurl_FCC( n1, n2, n3, Lambda_1, Lambda_2, Lambda_3 );
        cluster = '110';
        SVD_curl.cluster110 = SVDSingleCurl_FCC_cluster( n1, n2, n3, Lambda_1, Lambda_2, Lambda_3, cluster );
        cluster = '101';
        SVD_curl.cluster101 = SVDSingleCurl_FCC_cluster( n1, n2, n3, Lambda_1, Lambda_2, Lambda_3, cluster );
        cluster = '011';
        SVD_curl.cluster011 = SVDSingleCurl_FCC_cluster( n1, n2, n3, Lambda_1, Lambda_2, Lambda_3, cluster );
        
        % mtx_C = Sparse_mtx_C_FCC(n1, n2, n3, vec_k, period_vec1, period_vec2, period_vec3, delta);
        mtx_C.cluster000 = SVD_curl.cluster000.Lambda_Pr * (sparse(diag(SVD_curl.Sigma_r)) * SVD_curl.cluster000.Lambda_Qr');
        % mtx_C.cluster000(abs(mtx_C.cluster000) < 1e-13) = 0;
        mtx_C.cluster110 = SVD_curl.cluster110.Lambda_Pr * (sparse(diag(SVD_curl.Sigma_r)) * SVD_curl.cluster110.Lambda_Qr');
        % mtx_C.cluster110(abs(mtx_C.cluster110) < 1e-13) = 0;
        mtx_C.cluster101 = SVD_curl.cluster101.Lambda_Pr * (sparse(diag(SVD_curl.Sigma_r)) * SVD_curl.cluster101.Lambda_Qr');
        % mtx_C.cluster101(abs(mtx_C.cluster101) < 1e-13) = 0;
        mtx_C.cluster011 = SVD_curl.cluster011.Lambda_Pr * (sparse(diag(SVD_curl.Sigma_r)) * SVD_curl.cluster011.Lambda_Qr');
        
        mtx_T_prod_vec  = @(input_vec)mtx_T_prod_vec_FCC(input_vec, FFT_parameter);
        mtx_TH_prod_vec = @(input_vec)mtx_TH_prod_vec_FCC(input_vec, FFT_parameter);

        
end

fprintf(fid,'dim = %8.0f; n = %8.0f; \n', 24*n, n );
fprintf('dim = %8.0f; \n', 24*n );

fprintf(fid,'vec_k = [ %24.16e; %24.16e; %24.16e ]; \n', vec_k(1,1), vec_k(2,1), vec_k(3,1));
fprintf('icase = %2.0f, icount = %2.0f \n', icase, icount);
fprintf(fid,'icase = %2.0f; icount = %2.0f; \n', icase, icount);
fprintf(fid,'gamma = %24.16e; \n',chirality_parameter);
        

s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s);


i_idx     = zeros(16*n,1);
j_idx     = zeros(16*n,1);
val       = zeros(16*n,1);

i_idx(1:8*n,1) = 1:8*n;
j_idx(1:8*n,1) = 8*n+1:16*n;
val(1:2*n,1)   = 1i./SVD_curl.Sigma_r;
val(2*n+1:4*n,1)   = 1i./SVD_curl.Sigma_r;
val(4*n+1:6*n,1)   = 1i./SVD_curl.Sigma_r;
val(6*n+1:8*n,1)   = 1i./SVD_curl.Sigma_r;

i_idx(8*n+1:16*n,1) = 8*n+1:16*n;
j_idx(8*n+1:16*n,1) = 1:8*n;
val(8*n+1:10*n,1)   = -1i./SVD_curl.Sigma_r;
val(10*n+1:12*n,1)   = -1i./SVD_curl.Sigma_r;
val(12*n+1:14*n,1)   = -1i./SVD_curl.Sigma_r;
val(14*n+1:16*n,1)   = -1i./SVD_curl.Sigma_r;

GEP_mtx_B = sparse(i_idx, j_idx, val, 16*n, 16*n);

eigensolver = 'Krylov'; %'Krylov'; %'SIRA';

switch eigensolver
    
    case 'SIRA'
        
        LSopt.tol         = 1.0e-3;
        LSopt.infoConvBeh = 'yes';
        
        tic;  
        initial_V   = randn(16*n,1) + 1i * randn(16*n,1); % ones(dim,1);
        initial_V   = initial_V / norm(initial_V);
        %

        stop_tolerance = 1.0e-12; %1.0e8 * eps / scale; 
        target_type    = 'RGTR'; %'RGTC'; 

        no_restart     = 35; %35;
        mtxdim         = 16 * n;
        CorEq          = 'SIRA'; %'JD'; %'SIRA';
        LSinfo.solver  = 'bicgstabl';%'bicgstabl'; %'bicgstabl'; %'minres';
        LSinfo.precond = 'no'; %'yes';
        flag_LU = 0;
        

       [conv_ew, conv_ev] = GEP_AB_Herm_JDSIRA_Driver (@(vec)mtx_prod_vec_range_chiral(vec, SVD_curl, ...
           mtx_data, @(x)mtx_TH_prod_vec(x), @(x)mtx_T_prod_vec(x) ), GEP_mtx_B, mtxdim, no_restart, ...
            ncv, stop_tolerance, initial_V, sigma, fid, target_type, CorEq, LSinfo); % @(x)solve_Minv_b( x, LU_precond ));


    case 'Krylov'
        inner_iter_count  = 0;
        outer_iter_count  = 0;
        tol_LinearSys     = 1.0e-2 * eps / delta(1)^2; %5.0e2 * eps / delta(1)^2;
        %tol_LinearSys     = 5.0e-4 * eps / delta(1)^2;
        linear_system_tol = 1e-14; %1e-14;
        pcg_iter_number   = 10^4;
        
        scale       = 2 * sqrt(1/delta(1)^2 + 1/delta(2)^2 + 1/delta(3)^2);
        opts.tol    = 1.0e-12; %1.0e4 * eps / scale; %eps / scale;
        opts.maxit  = 1000;
        opts.p      = 25; %3 * nstp; %3 * nstp;  %5*nstp;
        opts.issym  = 0;
        opts.isreal = 0;
        opts.v0     = randn(24*n,1) + 1i * randn(24*n,1);
        
        [ev, ew] = eigs(@(vec)mtx_prod_vec_SEP(vec, mtx_C, mtx_data, n), 24*n, ncv, 'lm', opts);
%         [ev, ew] = eigs (@(vec)mtx_prod_vec_range_chiral(vec, SVD_curl, ...
%            mtx_data, @(x)mtx_TH_prod_vec(x), @(x)mtx_T_prod_vec(x)), 16*n, GEP_mtx_B, ncv, 'lr', opts);
%         res = mtx_prod_vec_SEP(ev(:,1), mtx_C, mtx_data, n) - ew(1) * ev(:,1);
%         norm(res);
     
        ew = ew / 1i;

end


fclose(fid);
