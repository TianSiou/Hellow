function [ NM_ele_left_matrix, NM_mag_left_matrix ] = FAME_Matrix_Boundary_Parameter_Nonuniform_Mesh( grid_number, grid_set, mesh_len )
% This routine return the diagonal entries of boundary matrices 'NM_ele_left_matrix' and 'NM_mag_left_matrix'
%   which transform a standard 1-d differential matrix into a nonuniform
%   mesh differential matrix.
% grid_number:
%   must be a positive integer which stand for the number of total
%   discretization area.
% grid_set: 
%   must be a array of length same as (grid_number+1) which stand positions 
%   for each grid.
% mesh_len:
%   must be a scalar which stand for length of uniform mesh.
    step_len_ele = ones(grid_num,1);
    step_len_mag = ones(grid_num,1);

    step_len_ele(1:grid_num)   = grid_set(2:grid_num) - grid_set(1:grid_num-1);
    step_len_mag(2:grid_num-1) = ( grid_set(3:grid_num) - grid_set(1:grid_num-2) )/2;
    step_len_mag(1)            = (grid_set(2) - grid_set(1) + grid_set(grid_num) - grid_set(grid_num-1))/2;
    step_len_mag(grid_num)     = (grid_set(2) - grid_set(1) + grid_set(grid_num) - grid_set(grid_num-1))/2;
    
    if grid_number == length(step_len);
        NM_ele_left_matrix = mesh_len./step_len_ele;
        NM_mag_left_matrix = mesh_len./step_len_mag;
    else
        warning('The length of step_len must equal to grid_number!')
    end

end