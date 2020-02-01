function Compute_EW_dep( wave_vec_idx )

warning off

fame_m_path_string = pwd;
fame_m_path_string = fame_m_path_string(1:end-15);
GTSHIRA_path_string = fullfile('My_routines','GTSHIRA');
addpath(fame_m_path_string)
addpath My_routines
addpath(GTSHIRA_path_string)

dir_name = fullfile(pwd,'result',[datestr(now,30),'k_',num2str(wave_vec_idx)]);
mkdir(dir_name);

FAME_Folder_Manager(fame_m_path_string) % loading path for whole FAME

%% Load user options from exist file "FAME_User_Option.m"
[ Popt ] = FAME_User_Option();

%% Generate modified lattice vectors and lattice constants for computing
[ Par.mesh, Par.lattice, Par.recip_lattice, Par.material, Par.eig ] = FAME_Parameter_Generator( Popt);
%% Generate wave vector array by partition the path of Brillouin zone 
path_string = 'PNHP';
[ Par.recip_lattice ] = FAME_Parameter_Brillouin_Zone_Path( Popt.recip_lattice.part_num, Par.lattice, Par.recip_lattice, path_string );
%% Locating indices for the material inside
[ Par.material.B.ele_x_idx, Par.material.B.ele_y_idx, Par.material.B.ele_z_idx, Par.material.B.mag_x_idx, Par.material.B.mag_y_idx, Par.material.B.mag_z_idx, Par.material.B.org_idx ] = ...
    FAME_Material_Locate_Index( Par.mesh, Par.lattice, Par.material);

%% Setting frequency and wave vector
wave_vec        = Par.recip_lattice.wave_vec_array(:,wave_vec_idx);
wave_vec        = wave_vec/norm(wave_vec);
freq_array      = linspace(0.35,0.6,100);
freq_array_2pi  = freq_array*2*pi;
%% Construct Lambdas and discrete curl-operators with periodic boundary conditions
[ C , Cs, C_1, C_2, C_3 ] = ...
    FAME_Matrix_Curl( [0;0;0], Par.mesh.grid_num, Par.mesh.edge_len,Par.mesh.mesh_len, {'quasi_periodic','quasi_periodic','quasi_periodic'}, Par.lattice.lattice_type, Par.lattice.lattice_vec_a, Par.lattice.lattice_constant );
[ Lambdas ] = ...
    FAME_Matrix_Lambdas( [0;0;0], Par.mesh.grid_num, Par.mesh.mesh_len, Par.lattice.lattice_type, Par.lattice.lattice_constant, Par.lattice.lattice_vec_a );
FAME_Check_Eigendecomp(Par.mesh.grid_num, Par.lattice.lattice_type, C, Cs, C_1, C_2, C_3, Lambdas);
%% Construct the wave vector dependent matrices K_tilde
N = Par.mesh.grid_num(1)*Par.mesh.grid_num(2)*Par.mesh.grid_num(3);
I_N  = speye(N);
O_N  = sparse(N, N);

K_1_tilde = .5*wave_vec(1)*(Par.mesh.mesh_len(1)*C_1 + 2*I_N);
K_2_tilde = .5*wave_vec(2)*(Par.mesh.mesh_len(2)*C_2 + 2*I_N);
K_3_tilde = .5*wave_vec(3)*(Par.mesh.mesh_len(3)*C_3 + 2*I_N);
K_tilde   = [        O_N, -K_3_tilde,  K_2_tilde;
               K_3_tilde,        O_N, -K_1_tilde;
              -K_2_tilde,  K_1_tilde,        O_N];

Lambdas.Lambda_x_tilde = .5*wave_vec(1)*(Par.mesh.mesh_len(1)*Lambdas.Lambda_x + 2);
Lambdas.Lambda_y_tilde = .5*wave_vec(2)*(Par.mesh.mesh_len(2)*Lambdas.Lambda_y + 2);
Lambdas.Lambda_z_tilde = .5*wave_vec(3)*(Par.mesh.mesh_len(3)*Lambdas.Lambda_z + 2);


for i = 1:length(freq_array_2pi) 
    omega_p = 9;
    gamma   = 0.054;
    Par.material.ele_permitt_in = 1 - (omega_p^2) / (freq_array(i)^2 - 1i*gamma*freq_array(i));
    
    B = FAME_Matrix_B_Isotropic( Par.mesh, Par.material);
    % Construst the coefficient matrices of GQEP, PQEP, and their transformations
    mtx = Matrix_Generate_1(C, C_1, C_2, C_3, K_tilde, K_1_tilde, K_2_tilde, K_3_tilde, B, Lambdas, N);
    
    fprintf('Start computing wave length with freqency(%d/%d) = %.3f,(2pi) and wave vector = [%.2f,%.2f,%.2f].\n',i,length(freq_array_2pi) ,freq_array(i),wave_vec);

    theta    = 0.5*pi;
    shift_nu0 = exp(1i*theta);
    fprintf('Use shift value nu_0 = exp(1i*%f*pi)\n',theta/pi);
    freq_2pi = freq_array_2pi(i);
    opt.linear_system = LS_settings( 3*N, shift_nu0, freq_array_2pi(i), B);
    opt.GTSHIRA       = GTSHIRA_settings( 3*N, shift_nu0 );
    
     %% Construst the coefficient matrices of GQEP
    % Construst the coefficient matrices of GQEP, PQEP, and their transformations
    mtx = Matrix_Generate_2(mtx, freq_array_2pi(i), shift_nu0, opt.linear_system.alpha );    
    % �غc�U�� function handle
    funhand = Funhand_Generate(C, K_tilde,B, mtx, Par.mesh, Par.lattice, Lambdas, freq_array_2pi(i), opt, N);
    % �p�� GQEP & PQEP ���S�x��
    [mu{i}, nu{i}, Ele_field, Info{i}] = compute_by_GTSHIRA(mtx, funhand, funhand.fun_invPp, opt);
    % �p�� GQEP & PQEP �S�x�諸�۹�~�t�P����~�t
    frobnorm_K = norm( -K_tilde.'*K_tilde, 'fro' );
    frobnorm_G = norm( C.'*K_tilde - K_tilde.'*C, 'fro' );
    frobnorm_M = norm( C.'*C - (freq_2pi^2)*spdiags(B.B_eps,0,3*N,3*N), 'fro' );
    frobnorm_A  = norm( 0.25*(C.'*C - (freq_2pi^2)*spdiags(B.B_eps,0,3*N,3*N) - (C.'*K_tilde - K_tilde.'*C) + (-K_tilde.'*K_tilde)), 'fro' );
    frobnorm_AT = norm( 0.25*(C.'*C - (freq_2pi^2)*spdiags(B.B_eps,0,3*N,3*N) + (C.'*K_tilde - K_tilde.'*C) + (-K_tilde.'*K_tilde)), 'fro' );
    frobnorm_Q  = norm( 0.5 *(C.'*C - (freq_2pi^2)*spdiags(B.B_eps,0,3*N,3*N) - (-K_tilde.'*K_tilde)), 'fro' );
    temp_g = @(lambda) abs(lambda)^2*frobnorm_M  + abs(lambda)*frobnorm_G + frobnorm_K;
    temp_p = @(lambda) abs(lambda)^2*frobnorm_AT - abs(lambda)*frobnorm_Q + frobnorm_A;
    res = {};
    for j = 1:length(mu{i})
        res.absres_p(j) = norm(funhand.palindromic.Qp(Ele_field(:,j),nu{i}(j)));
        res.relres_p(j) = res.absres_p(j)/(temp_p(nu{i}(j))*norm(Ele_field(:,j)));
        res.absres_g(j) = norm(funhand.palindromic.Qg(Ele_field(:,j),mu{i}(j)));
        res.relres_g(j) = res.absres_g(j)/(temp_g(mu{i}(j))*norm(Ele_field(:,j)));
    end
    lambda = 1i./mu{i};
    figure(1); plot(freq_array(i)*ones(length(lambda),1),real(lambda),'r*'); hold on
    figure(2); plot(freq_array(i)*ones(length(lambda),1),imag(lambda),'b*'); hold on
    
    file_name = fullfile(dir_name,['result_k',num2str(wave_vec_idx),'_w',num2str(i)]);
    Info{i}.linear_system.v0 = [];
    Info{i}.GTSHIRA.v0       = [];
    save(file_name,'mu','nu','res','Info');
    fprintf('Done! save results as file ')
end
end
function plot_unit_circle(a,b,c)
theta = 0:2*pi/100:2*pi;
x = cos(theta);
y = sin(theta);

subplot(a,b,c); hold on;
plot(x,y,'k-');    
hold off;
end