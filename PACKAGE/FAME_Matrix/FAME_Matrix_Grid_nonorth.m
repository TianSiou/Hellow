function [ Origin_point_set ] = FAME_Matrix_Grid_nonorth( grid_nums, lattice_vec_a )
    n_grid_num = grid_nums(1)*grid_nums(2)*grid_nums(3);
    
    %% Create the Yee's scheme point
    % Creating the grids on the vertices
    x_grid  = (   0 : grid_nums(1)-  1 )' / grid_nums(1);
    y_grid  = (   0 : grid_nums(2)-  1 )' / grid_nums(2);
    z_grid  = (   0 : grid_nums(3)-  1 )' / grid_nums(3);
    %% Constructing points on edge for the electric field and on face for the magnetic field
%     Origin_point_set     = zeros(n_grid_num , 3);
    % Constructing points on original grid
    Origin_point_set = [ kron( ones(grid_nums(3)*grid_nums(2) , 1) ,                                    x_grid  ) ,...
                         kron( ones(grid_nums(3) ,         1) , kron(       y_grid , ones(grid_nums(1) , 1)  )  ) ,...
                         kron(                         z_grid ,             ones(grid_nums(2)*grid_nums(1) , 1) )  ] * lattice_vec_a';
end