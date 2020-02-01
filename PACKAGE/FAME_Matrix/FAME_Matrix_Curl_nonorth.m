function [ C , Cs, C_1, C_2, C_3, C_x, C_y, C_z ] = FAME_Matrix_Curl_nonorth( wave_vec, grid_nums, edge_lens, mesh_lens, bd_conds, lattice_type, lattice_vec_a, lattice_constant, pml_grid_numbers )
% this routine return the matrices of discrete curl operators C and Cs
% wave_vec: 
%   a desired 3-d vector.
% grid_nums: 
%   a 1x3 matrix which contained the grid numbers of x,y,z directions in the
%   discretization.
% bd_conds:
%   a cell which contained 3 strings which stand for these desired boundary 
%   conditions in x,y,z directions. These string can be given between 
%                'quasi_periodic', 'dirichlet' and 'perfect_matching_layer'
% lattice_type:
%   a string which stand for the desired lattice type. This string can be 
%   given between 'simple_cubic', 'face_centered_cubic' and 'body_centered_cubic'
% lattice_constant:
%   a positive number which stand for the desired lattice constant.
% layer_thickness:(must be given as bd_conds has 'perfect_matching_layer')
%   a 1x3 matrix which obtained the desired thickness of each perfect matching
%   layer.
% pml_grid_numbers
%   a 1x3 matrix which contained the grid numbers of x,y,z directions in the
%   each perfect matching layer.


idouble_pi = 1i*2*pi;
n_grid_num = grid_nums(1)*grid_nums(2)*grid_nums(3);
one_x      = ones(grid_nums(1),1);
one_y      = ones(grid_nums(2),1);
one_z      = ones(grid_nums(3),1);

lattice_vec_a1 = lattice_vec_a(:,1);
lattice_vec_a2 = lattice_vec_a(:,2);
lattice_vec_a3 = lattice_vec_a(:,3);

theta_x = dot(lattice_vec_a1,wave_vec);
theta_y = dot(lattice_vec_a2,wave_vec);
theta_z = dot(lattice_vec_a3,wave_vec);

switch lattice_type
    case {'simple_cubic', 'primitive_orthorhombic', 'primitive_tetragonal'}
        I_x   = speye(grid_nums(1),grid_nums(1));
        I_y   = speye(grid_nums(2),grid_nums(2));
        I_z   = speye(grid_nums(3),grid_nums(3));
        %% Constructing the differential matrices and bounday matrices for different boundary conditions
        K_x   = spdiags([-one_x one_x],0:1,grid_nums(1),grid_nums(1));
        K_y   = spdiags([-one_y one_y],0:1,grid_nums(2),grid_nums(2));
        K_z   = spdiags([-one_z one_z],0:1,grid_nums(3),grid_nums(3));

        % Correcting the differential matrices for different boundary conditions
        switch bd_conds{1}    
            case 'dirichlet'
                  K_xs = K_x';
            case 'quasi_periodic'
                  K_x(end,1) = exp(idouble_pi*theta_x);
                  K_xs = K_x';
            case 'perfect_matching_layer'
%                 [ pml_ele_left_matrix_x, pml_mag_left_matrix_x ] = ...
%                     FAME_Construct_Matrix_Boundary_Parameter_PML( grid_nums(1), pml_grid_numbers(1), mesh_lens(1) );
%                 K_x_temp = K_x;
%                 K_x  = spdiags(pml_ele_left_matrix_x,0,grid_nums(1),grid_nums(1))*K_x_temp ;
%                 K_xs = spdiags(pml_mag_left_matrix_x,0,grid_nums(1),grid_nums(1))*K_x_temp';
                pml.x = FAME_Construct_Matrix_Boundary_Parameter_PML( grid_nums(1), pml_grid_numbers(1), mesh_lens(1), edge_lens(1) );
                K_x_temp = K_x;
                K_x  = spdiags(pml_ele_left_matrix_x,0,grid_nums(1),grid_nums(1))*K_x_temp ;
                K_xs = spdiags(pml_mag_left_matrix_x,0,grid_nums(1),grid_nums(1))*K_x_temp';
        end
        switch bd_conds{2}    
            case 'dirichlet'
                  K_ys = K_y';
            case 'quasi_periodic'
                  K_y(end,1) = exp(idouble_pi*theta_y);
                  K_ys = K_y';
            case 'perfect_matching_layer'
                [ pml_ele_left_matrix_y, pml_mag_left_matrix_y ] = ...
                    FAME_Construct_Matrix_Boundary_Parameter_PML( grid_nums(2), pml_grid_numbers(2), mesh_lens(2) );
                K_y_temp = K_y;
                K_y  = spdiags(pml_ele_left_matrix_y,0,grid_nums(2),grid_nums(2))*K_y_temp ;
                K_ys = spdiags(pml_mag_left_matrix_y,0,grid_nums(2),grid_nums(2))*K_y_temp';
        end
        switch bd_conds{3}    
            case 'dirichlet'
                  K_zs = K_z';
            case 'quasi_periodic'
                  K_z(end,1) = exp(idouble_pi*theta_z);
                  K_zs = K_z';
            case 'perfect_matching_layer'
                [ pml_ele_left_matrix_z, pml_mag_left_matrix_z ] = ...
                    FAME_Construct_Matrix_Boundary_Parameter_PML( grid_nums(3), pml_grid_numbers(3), mesh_lens(3) );
                K_z_temp = K_z;
                K_z  = spdiags(pml_ele_left_matrix_z,0,grid_nums(3),grid_nums(3))*K_z_temp ;
                K_zs = spdiags(pml_mag_left_matrix_z,0,grid_nums(3),grid_nums(3))*K_z_temp';
        end
        C_1  =  kron( I_z , kron( I_y , K_x ) )./mesh_lens(1);
        C_2  =  kron( I_z , kron( K_y , I_x ) )./mesh_lens(2);
        C_3  =  kron( K_z , kron( I_y , I_x ) )./mesh_lens(3);
        Ce_1 = C_1 + ...
               kron( kron( I_z, K_y - K_y'), .5*(K_x + 2*I_x)) / mesh_lens(2) + ...
               kron( kron( K_z - K_z', I_y), .5*(K_x + 2*I_x)) / mesh_lens(3);
        Ce_2 = kron( kron( I_z, .5*(K_y + 2*I_y)), K_x - K_x') / mesh_lens(1) + ...
               C_2 + ...
               kron( kron( K_z - K_z', .5*(K_y + 2*I_y)), I_x) / mesh_lens(3);
        Ce_3 = kron( kron( .5*(K_z + 2*I_z), I_y), K_x - K_x') / mesh_lens(1) + ...
               kron( kron( .5*(K_z + 2*I_z), K_y - K_y'), I_x) / mesh_lens(2) + ...
               C_3;
        Va = det(lattice_vec_a);   
        Ce_1 = Va*Ce_1/(2*pi);
        Ce_2 = Va*Ce_2/(2*pi);
        Ce_3 = Va*Ce_3/(2*pi);

        C_1s =  kron( I_z  , kron( I_y  , K_xs ) )./mesh_lens(1);
        C_2s =  kron( I_z  , kron( K_ys , I_x  ) )./mesh_lens(2);
        C_3s =  kron( K_zs , kron( I_y  , I_x  ) )./mesh_lens(3);
        Ce_1s = Ce_1';
        Ce_2s = Ce_2';
        Ce_3s = Ce_3';
    otherwise
        n1 = grid_nums(1);         n2 = grid_nums(2);         n3 = grid_nums(3);
        m1 = lattice_constant.m1;  m2 = lattice_constant.m2;  m3 = lattice_constant.m3;  m4 = lattice_constant.m4;

        i2pika1   =  idouble_pi*dot(wave_vec,lattice_vec_a1);
        i2pika2   =  idouble_pi*dot(wave_vec,lattice_vec_a2);
        i2pika3   =  idouble_pi*dot(wave_vec,lattice_vec_a3);

        I_x    = speye(grid_nums(1));
        I_y    = speye(grid_nums(2));
        I_xy   = speye(grid_nums(1)*grid_nums(2));
        I_z    = speye(grid_nums(3));
        I_m3   = speye(m3);
        I_n2m3 = speye(n2-m3);
        
        K_1   =  spdiags([-one_x one_x],0:1,grid_nums(1),grid_nums(1));

        k_2    = spdiags([-one_y one_y],0:1,grid_nums(2),grid_nums(2));
        K_2    = kron(k_2,I_x);

        k_3    = spdiags([-one_z one_z],0:1,grid_nums(3),grid_nums(3));
        K_3    = kron(k_3,I_xy);

        switch bd_conds{1}  
            case 'dirichlet'
                  K_1s = K_1';
            case 'quasi_periodic'
                  K_1(end,1) = exp(i2pika1);
                  K_1s = K_1';
            case 'perfect_matching_layer'
%                 [ pml_ele_left_matrix_x, pml_mag_left_matrix_x ] = ...
%                     FAME_Construct_Matrix_Boundary_Parameter_PML( grid_nums(1), pml_grid_numbers(1), mesh_lens(1) );
%                 K_1_temp = K_1;
%                 K_1  = spdiags(pml_ele_left_matrix_x,0,grid_nums(1),grid_nums(1))*K_1_temp ;
%                 K_1s = spdiags(pml_mag_left_matrix_x,0,grid_nums(1),grid_nums(1))*K_1_temp';
        end
        switch bd_conds{2} 
            case 'dirichlet'
                  K_2s = K_2';
            case 'quasi_periodic'
                  if lattice_constant.theta_3 <= pi/2
                      rho_1 = 0;
                  else
                      rho_1 = 1;
                  end
                  idx_i     = 1:n1;
                  idx_j     = [n1-m1+1:n1  1:n1-m1];
                  s         = exp(i2pika1*rho_1)*[exp(-i2pika1)*ones(1,m1) ones(1,n1-m1)];
                  J_2       = sparse(idx_i,idx_j,s,n1,n1);
                  
                  K_2(n1*(n2-1)+1:n1*n2,1:n1) = exp(i2pika2)*J_2;

                  K_2s = K_2';
            case 'perfect_matching_layer'
%                 [ pml_ele_left_matrix_y, pml_mag_left_matrix_y ] = ...
%                     FAME_Construct_Matrix_Boundary_Parameter_PML( grid_nums(2), pml_grid_numbers(2), mesh_lens(2) );
%                 K_2_temp = K_2;
%                 K_2  = spdiags( kron(pml_ele_left_matrix_y,ones(grid_nums(1),1)), 0,grid_nums(2)*grid_nums(1),grid_nums(2)*grid_nums(1))*K_2_temp ;
%                 K_2s = spdiags( kron(pml_mag_left_matrix_y,ones(grid_nums(1),1)), 0,grid_nums(2)*grid_nums(1),grid_nums(2)*grid_nums(1))*K_2_temp';
        end
        switch bd_conds{3}  
            case 'dirichlet'
                  K_3s = K_3';
            case 'quasi_periodic'
                  i2pikt1 = idouble_pi*dot(wave_vec,lattice_constant.t1);
                  i2pikt2 = idouble_pi*dot(wave_vec,lattice_constant.t2);
                  i2pikt3 = idouble_pi*dot(wave_vec,lattice_constant.t3);
                  i2pikt4 = idouble_pi*dot(wave_vec,lattice_constant.t4);

                  idx_i     = 1:n1;
                  idx_j     = [n1-m4+1:n1  1:n1-m4];
                  s         = [exp(i2pikt3)*ones(1,m4) exp(i2pikt4)*ones(1,n1-m4)];
                  temp1     = kron(I_m3, sparse(idx_i,idx_j,s,n1,n1) );
                  
                  idx_i     = 1:n1;
                  idx_j     = [n1-m2+1:n1  1:n1-m2];
                  s         = [exp(i2pikt2)*ones(1,m2) exp(i2pikt1)*ones(1,n1-m2)];
                  temp2     = kron(I_n2m3, sparse(idx_i,idx_j,s,n1,n1) );
% exp(i2pikt1)
% exp(i2pikt2)
% exp(i2pikt3)
% exp(i2pikt4)
% t1 = lattice_constant.t1
% t2 = lattice_constant.t2
% t3 = lattice_constant.t3
% t4 = lattice_constant.t4
                  O_1 = sparse( n1*m3, n1*(n2-m3) );
                  O_2 = sparse( n1*(n2-m3), n1*m3 );
                  J_3 = [O_1, temp1; temp2, O_2];
% figure(2); cla
% spy(J_3)
                  K_3(n1*n2*(n3-1)+1:n1*n2*n3,1:n1*n2) = exp(i2pika3)*J_3;
                  K_3s = K_3';
            case 'perfect_matching_layer'
%                 [ pml_ele_left_matrix_z, pml_mag_left_matrix_z ] = ...
%                     FAME_Construct_Matrix_Boundary_Parameter_PML( grid_nums(3), pml_grid_numbers(3), mesh_lens(3) );
%                 K_3_temp = K_3;
%                 K_3  = spdiags( kron(pml_ele_left_matrix_z,ones(grid_nums(2)*grid_nums(1),1)), 0,grid_nums(3)*grid_nums(2)*grid_nums(1),grid_nums(3)*grid_nums(2)*grid_nums(1))*K_3_temp ;
%                 K_3s = spdiags( kron(pml_mag_left_matrix_z,ones(grid_nums(2)*grid_nums(1),1)), 0,grid_nums(3)*grid_nums(2)*grid_nums(1),grid_nums(3)*grid_nums(2)*grid_nums(1))*K_3_temp';
        end

        C_1  =  kron( I_z , kron( I_y , K_1 ) )/mesh_lens(1);
        C_2  =  kron( I_z , K_2 )/mesh_lens(2);
        C_3  =  K_3/mesh_lens(3);

        C_1s  =  kron( I_z , kron( I_y , K_1s ) )/mesh_lens(1);
        C_2s  =  kron( I_z , K_2s )/mesh_lens(2);
        C_3s  =  K_3s/mesh_lens(3);    
end
%% Constructing the curl matrix C and C*
  % Constructing the curl matrix C
O = sparse(n_grid_num,n_grid_num);
C = [    O, -C_3,  C_2;
       C_3,    O, -C_1;
      -C_2,  C_1,    O];
Cs = C';  
I  = speye(n_grid_num);
C_x = (1/6)*((mesh_lens(1)*C_1 + 2*I)*(mesh_lens(1)*C_1s + 2*I) + 2*I);
C_y = (1/6)*((mesh_lens(2)*C_2 + 2*I)*(mesh_lens(2)*C_2s + 2*I) + 2*I);
C_z = (1/6)*((mesh_lens(3)*C_3 + 2*I)*(mesh_lens(3)*C_3s + 2*I) + 2*I);
% eig(full((C_x)))