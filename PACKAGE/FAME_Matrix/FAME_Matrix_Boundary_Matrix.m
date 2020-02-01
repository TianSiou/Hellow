function [ pml, numesh ] = FAME_Matrix_Boundary_Matrix( grid_nums, pml_grid_nums, mesh_lens, edge_lens )
    N     = grid_nums(1)*grid_nums(2)*grid_nums(3);
    one_x = ones(grid_nums(1),1);
    one_y = ones(grid_nums(2),1);
    one_z = ones(grid_nums(3),1);
    eye_x = speye(grid_nums(1));
    eye_y = speye(grid_nums(2));
    eye_z = speye(grid_nums(3));

    %% Construct 1D boundary transform matrices in x,y,z-terms
    [ pml.x, numesh.x ] = FAME_Construct_Matrix_Boundary_Parameter_PML( grid_nums(1), pml_grid_nums(1), mesh_lens(1), edge_lens(1) );
    [ pml.y, numesh.y ] = FAME_Construct_Matrix_Boundary_Parameter_PML( grid_nums(2), pml_grid_nums(2), mesh_lens(2), edge_lens(2) );
    [ pml.z, numesh.z ] = FAME_Construct_Matrix_Boundary_Parameter_PML( grid_nums(3), pml_grid_nums(3), mesh_lens(3), edge_lens(3) );
    %% Construct 3D boundary transform matrices
    % x
    pml.ele_left_matrix_1      = kron( one_z, kron( one_y, pml.x.ele_left_matrix      ) );
    pml.ele_left_matrix_1_inv  = kron( one_z, kron( one_y, pml.x.ele_left_matrix_inv  ) );
    pml.mag_left_matrix_1      = kron( one_z, kron( one_y, pml.x.mag_left_matrix      ) );
    pml.mag_left_matrix_1_inv  = kron( one_z, kron( one_y, pml.x.mag_left_matrix_inv  ) );
%     pml.ele_left_matrix_1      = spdiags( kron( one_z, kron( one_y, pml.x.ele_left_matrix      ) ), 0, N, N )...
%                                               - kron(      eye_z, kron(      eye_y, pml.x.ele_left_lowrank     ) );
%     pml.ele_left_matrix_1_inv  = spdiags( kron( one_z, kron( one_y, pml.x.ele_left_matrix_inv  ) ), 0, N, N )...
%                                               + kron(      eye_z, kron(      eye_y, pml.x.ele_left_lowrank_inv ) );
%     pml.mag_left_matrix_1      = spdiags( kron( one_z, kron( one_y, pml.x.mag_left_matrix      ) ), 0, N, N );
%     pml.mag_left_matrix_1_inv  = spdiags( kron( one_z, kron( one_y, pml.x.mag_left_matrix_inv  ) ), 0, N, N );
% 
%     pml.mag_right_matrix_1     =          kron(      eye_z, kron(      eye_y, pml.x.mag_right_matrix     ) );
%     pml.mag_right_matrix_1_inv =          kron(      eye_z, kron(      eye_y, pml.x.mag_right_matrix_inv ) );
    % y
    pml.ele_left_matrix_2      = kron( one_z, kron( pml.y.ele_left_matrix     , one_x ) );
    pml.ele_left_matrix_2_inv  = kron( one_z, kron( pml.y.ele_left_matrix_inv , one_x ) );
    pml.mag_left_matrix_2      = kron( one_z, kron( pml.y.mag_left_matrix     , one_x ) );
    pml.mag_left_matrix_2_inv  = kron( one_z, kron( pml.y.mag_left_matrix_inv , one_x ) );
%     pml.ele_left_matrix_2      = spdiags( kron( one_z, kron( pml.y.ele_left_matrix     , one_x ) ), 0, N, N )...
%                                               - kron(      eye_z, kron( pml.y.ele_left_lowrank    ,      eye_x ) );
%     pml.ele_left_matrix_2_inv  = spdiags( kron( one_z, kron( pml.y.ele_left_matrix_inv , one_x ) ), 0, N, N )...
%                                                   + kron(      eye_z, kron( pml.y.ele_left_lowrank_inv,      eye_x ) );
%     pml.mag_left_matrix_2      = spdiags( kron( one_z, kron( pml.y.mag_left_matrix     , one_x ) ), 0, N, N );
%     pml.mag_left_matrix_2_inv  = spdiags( kron( one_z, kron( pml.y.mag_left_matrix_inv , one_x ) ), 0, N, N );
% 
%     pml.mag_right_matrix_2     =          kron(      eye_z, kron( pml.y.mag_right_matrix    ,      eye_x ) );
%     pml.mag_right_matrix_2_inv =          kron(      eye_z, kron( pml.y.mag_right_matrix_inv,      eye_x ) );
    % z
    pml.ele_left_matrix_3      = kron( pml.z.ele_left_matrix     , kron( one_y, one_x ) );
    pml.ele_left_matrix_3_inv  = kron( pml.z.ele_left_matrix_inv , kron( one_y, one_x ) );
    pml.mag_left_matrix_3      = kron( pml.z.mag_left_matrix     , kron( one_y, one_x ) );
    pml.mag_left_matrix_3_inv  = kron( pml.z.mag_left_matrix_inv , kron( one_y, one_x ) );
%     pml.ele_left_matrix_3      = spdiags( kron( pml.z.ele_left_matrix     , kron( one_y, one_x ) ), 0, N, N )...
%                                               - kron( pml.z.ele_left_lowrank    , kron(      eye_y,      eye_x ) );
%     pml.ele_left_matrix_3_inv  = spdiags( kron( pml.z.ele_left_matrix_inv , kron( one_y, one_x ) ), 0, N, N )...
%                                                   + kron( pml.z.ele_left_lowrank_inv, kron(      eye_y,      eye_x ) );
%     pml.mag_left_matrix_3      = spdiags( kron( pml.z.mag_left_matrix     , kron( one_y, one_x ) ), 0, N, N );
%     pml.mag_left_matrix_3_inv  = spdiags( kron( pml.z.mag_left_matrix_inv , kron( one_y, one_x ) ), 0, N, N );
% 
%     pml.mag_right_matrix_3     =          kron( pml.z.mag_right_matrix    , kron(      eye_y,      eye_x ) );
%     pml.mag_right_matrix_3_inv =          kron( pml.z.mag_right_matrix_inv, kron(      eye_y,      eye_x ) );
pml.L_ele = [ pml.ele_left_matrix_2 .* pml.ele_left_matrix_3;...
              pml.ele_left_matrix_1 .* pml.ele_left_matrix_3;...
              pml.ele_left_matrix_1 .* pml.ele_left_matrix_2 ];
pml.L_mag = [ pml.mag_left_matrix_2 .* pml.mag_left_matrix_3;...
              pml.mag_left_matrix_1 .* pml.mag_left_matrix_3;...
              pml.mag_left_matrix_1 .* pml.mag_left_matrix_2 ];
pml.R_ele = [ pml.ele_left_matrix_1_inv; pml.ele_left_matrix_2_inv; pml.ele_left_matrix_3_inv ];
pml.R_mag = [ pml.mag_left_matrix_1_inv; pml.mag_left_matrix_2_inv; pml.mag_left_matrix_3_inv ];

pml.L_ele_inv = [ pml.ele_left_matrix_2_inv .* pml.ele_left_matrix_3_inv;...
                  pml.ele_left_matrix_1_inv .* pml.ele_left_matrix_3_inv;...
                  pml.ele_left_matrix_1_inv .* pml.ele_left_matrix_2_inv ];
pml.L_mag_inv = [ pml.mag_left_matrix_2_inv .* pml.mag_left_matrix_3_inv;...
                  pml.mag_left_matrix_1_inv .* pml.mag_left_matrix_3_inv;...
                  pml.mag_left_matrix_1_inv .* pml.mag_left_matrix_2_inv ];
pml.R_ele_inv = [ pml.ele_left_matrix_1; pml.ele_left_matrix_2; pml.ele_left_matrix_3 ];
pml.R_mag_inv = [ pml.mag_left_matrix_1; pml.mag_left_matrix_2; pml.mag_left_matrix_3 ];

end