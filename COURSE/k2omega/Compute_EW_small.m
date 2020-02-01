function Compute_EW( wave_vec_idx )

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

B = FAME_Matrix_B_Isotropic( Par.mesh, Par.material);

% Construst the coefficient matrices of GQEP, PQEP, and their transformations
% mtx = Matrix_Generate_1(C, C_1, C_2, C_3, K_tilde, K_1_tilde, K_2_tilde, K_3_tilde, B, Lambdas, N);
% for i = 1:length(freq_array_2pi) 
for i = 1:1
    mtx_gyroscopic  = Matrix_Gyroscopic( C, K_tilde, B );    
    mtx_gyroscopic.M = mtx_gyroscopic.fun_M(freq_array(i));
    % 在小矩陣時畫出各種轉化下的 Spectrum
    EW_GQEP  = polyeig( mtx_gyroscopic.K, mtx_gyroscopic.G, mtx_gyroscopic.M);
    EW_PQEP  = (1+EW_GQEP)./(1-EW_GQEP);
    EW_SS    = EW_PQEP + 1./EW_PQEP;
    lambda   = -1i*(1./EW_GQEP);
    
    idx_1 = find(abs(lambda)>1e+3);
    idx_2 = find(abs(lambda)<1e+3);
    idx_3 = intersect( find(abs(imag(lambda))<1e-2),idx_2);
    
    fprintf('Done!\n');
    figure(1); hold on
    plot(lambda(idx_1),'k.','markersize',7);
    plot(lambda(idx_2),'b*');
    axis equal
    figure(2); hold on
    plot(lambda(idx_2),'b*');
    plot(lambda(idx_3),'ro');
    axis([-40,40,-55,55])
    figure(3); hold on
    plot(EW_GQEP(idx_1),'k.');
    plot(EW_GQEP(idx_2),'b*');
    plot(EW_GQEP(idx_3),'ro');
    axis equal
    figure(4); hold on
    plot(EW_PQEP(idx_1),'k.');
    plot(EW_PQEP(idx_2),'b*');
    plot(EW_PQEP(idx_3),'ro');
    plot_unit_circle(1,4,3);
    axis equal
    figure(5); hold on
    plot(EW_SS(idx_1),'k.');
    plot(EW_SS(idx_2),'b*');
    plot(EW_SS(idx_3),'ro');
    plot_unit_circle(1,4,4);
    axis equal
    
    pause
end
end
function plot_unit_circle(a,b,c)
theta = 0:2*pi/100:2*pi;
x = cos(theta);
y = sin(theta);

% subplot(a,b,c); hold on;
plot(x,y,'k-');    
hold off;
end