function [ optimal_factor, optimal_grid, optimal_mesh_length, optimal_thickness, optimal_thickness_ratio ] = FAME_Boundary_Perfect_Matching_Layer_Optimal_Mesh( grid_num, layer_grid_num, layer_thickness, mesh_len, pml_opt )
max_loss_parameter = 8*( pml_opt.power+1 ) / ( pml_opt.frequency*layer_thickness );
optimal_factor     = 0.5*sqrt( 1+max_loss_parameter*max_loss_parameter )*grid_num/layer_grid_num;
sum_mesh_length    = layer_thickness;

optimal_factor_left  = 0;
optimal_factor_right = sqrt( 1+max_loss_parameter*max_loss_parameter )*grid_num/layer_grid_num;

while abs(sum_mesh_length) > 1e-12
optimal_thickness_ratio = zeros(layer_grid_num, 1);
optimal_mesh_length     = zeros(layer_grid_num, 1);
sum_mesh_length         = layer_thickness;

for grid_idx = 1:layer_grid_num
optimal_thickness_ratio(grid_idx) = sum_mesh_length/layer_thickness;
optimal_mesh_length(grid_idx)     = optimal_factor * mesh_len/sqrt( 1+max_loss_parameter*max_loss_parameter*( optimal_thickness_ratio(grid_idx)^(2*pml_opt.power) ) );
sum_mesh_length                   = sum_mesh_length - optimal_mesh_length(grid_idx);
end

if sum_mesh_length > 0
    optimal_factor_left = optimal_factor;
    optimal_factor = 0.5*(optimal_factor + optimal_factor_right);
else
    optimal_factor_right = optimal_factor;
    optimal_factor = 0.5*(optimal_factor_left + optimal_factor);
end

% if sum_mesh_length > 0
% optimal_factor = optimal_factor + sum_mesh_length;
% else
% optimal_factor = optimal_factor - sum_mesh_length;
% end
end

optimal_grid = cumsum(optimal_mesh_length);
optimal_grid = [0; optimal_grid];
optimal_thickness = optimal_thickness_ratio*layer_thickness;
optimal_factor = 1/optimal_factor;
end