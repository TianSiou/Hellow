function [ Electric_x_point_set, Electric_y_point_set, Electric_z_point_set,...
           Magnetic_x_point_set, Magnetic_y_point_set, Magnetic_z_point_set,...
               Origin_point_set ] = FAME_Matrix_Grid( grid_nums, mesh_lens )
    n_grid_num = grid_nums(1)*grid_nums(2)*grid_nums(3);
    
    %% Create the Yee's scheme point
    % Creating the grids on the vertices
    x_grid  = (   0 : grid_nums(1)-  1 ).' *mesh_lens(1);
    y_grid  = (   0 : grid_nums(2)-  1 ).' *mesh_lens(2);
    z_grid  = (   0 : grid_nums(3)-  1 ).' *mesh_lens(3);
    % Creating the grids on the edge center
    x_shift_grid = ( 0.5 : grid_nums(1)-0.5 ).' *mesh_lens(1);
    y_shift_grid = ( 0.5 : grid_nums(2)-0.5 ).' *mesh_lens(2);
    z_shift_grid = ( 0.5 : grid_nums(3)-0.5 ).' *mesh_lens(3);
    %% Constructing points on edge for the electric field and on face for the magnetic field
    Origin_point_set   = zeros(n_grid_num , 3);
    Electric_x_point_set = zeros(n_grid_num , 3);
    Electric_y_point_set = zeros(n_grid_num , 3);
    Electric_z_point_set = zeros(n_grid_num , 3);
    Magnetic_x_point_set = zeros(n_grid_num , 3);
    Magnetic_y_point_set = zeros(n_grid_num , 3);
    Magnetic_z_point_set = zeros(n_grid_num , 3);
    % Constructing points on original grid
    Origin_point_set(:,1) = kron( ones(grid_nums(3)*grid_nums(2) , 1) ,       x_grid );
    Origin_point_set(:,2) = kron( ones(grid_nums(3) , 1) , kron(       y_grid , ones(grid_nums(1) , 1)  )  );
    Origin_point_set(:,3) = kron(       z_grid , ones(grid_nums(2)*grid_nums(1) , 1) );
    % Constructing points on edge for the electric field
    Electric_x_point_set(:,1) = kron( ones(grid_nums(3)*grid_nums(2) , 1) , x_shift_grid );
    Electric_y_point_set(:,1) = kron( ones(grid_nums(3)*grid_nums(2) , 1) ,       x_grid );
    Electric_z_point_set(:,1) = kron( ones(grid_nums(3)*grid_nums(2) , 1) ,       x_grid );

    Electric_x_point_set(:,2) = kron( ones(grid_nums(3) , 1) , kron(       y_grid , ones(grid_nums(1) , 1)  )  );
    Electric_y_point_set(:,2) = kron( ones(grid_nums(3) , 1) , kron( y_shift_grid , ones(grid_nums(1) , 1)  )  );
    Electric_z_point_set(:,2) = kron( ones(grid_nums(3) , 1) , kron(       y_grid , ones(grid_nums(1) , 1)  )  ); 

    Electric_x_point_set(:,3) = kron(       z_grid , ones(grid_nums(2)*grid_nums(1) , 1) );
    Electric_y_point_set(:,3) = kron(       z_grid , ones(grid_nums(2)*grid_nums(1) , 1) );
    Electric_z_point_set(:,3) = kron( z_shift_grid , ones(grid_nums(2)*grid_nums(1) , 1) );
    % Constructing points on face for the magnetic field
    Magnetic_x_point_set(:,1) = kron( ones(grid_nums(3)*grid_nums(2) , 1) ,      x_grid );
    Magnetic_y_point_set(:,1) = kron( ones(grid_nums(3)*grid_nums(2) , 1) ,x_shift_grid );
    Magnetic_z_point_set(:,1) = kron( ones(grid_nums(3)*grid_nums(2) , 1) ,x_shift_grid );

    Magnetic_x_point_set(:,2) = kron( ones(grid_nums(3) , 1) , kron(y_shift_grid , ones(grid_nums(1) , 1)  )  );
    Magnetic_y_point_set(:,2) = kron( ones(grid_nums(3) , 1) , kron(      y_grid , ones(grid_nums(1) , 1)  )  );
    Magnetic_z_point_set(:,2) = kron( ones(grid_nums(3) , 1) , kron(y_shift_grid , ones(grid_nums(1) , 1)  )  ); 

    Magnetic_x_point_set(:,3) = kron(z_shift_grid , ones(grid_nums(2)*grid_nums(1) , 1) );
    Magnetic_y_point_set(:,3) = kron(z_shift_grid , ones(grid_nums(2)*grid_nums(1) , 1) );
    Magnetic_z_point_set(:,3) = kron(      z_grid , ones(grid_nums(2)*grid_nums(1) , 1) );
end