function  [ pml, numesh ] = FAME_Matrix_Boundary_Parameter_PML( grid_num, pml_grid_num, mesh_len, edge_len )
% This routine return the diagonal entries of boundary matrices 'pml_ele_left_matrix' and 'pml_mag_left_matrix'
%   which transform a standard 1-d differential matrix into a pml
%   differential matrix.
% grid_number:
%   must be a positive integer which stand for the number of total
%   discretization area.
% pml_grid_number: 
%   must be a positive integer which stand for the number of grids in each PML.
% inner_mesh_len:
%   must be a scalar which stand for length of mesh between the PMLs.
    pml_opt.power     = 4;
    pml_opt.ratio     = 1e-2;
    pml_opt.frequency = 0.5;
    layer_thickness   = pml_grid_num*mesh_len;
    inner_grid_num    = grid_num - 2*pml_grid_num;
    if inner_grid_num >= 0
        % Setting the inner value for Perfect Matching Layer matrix
        inner_one = ones(inner_grid_num,1);
            %% Mesh Optimization for PML Layers
        % ¢u¡C¢q¡C¢q¡C¢q¡C¢q¡@¢q¡@¢q¡@¢q¡@¢q¡C¢q¡C¢q¡C¢q¡C¢t
        % ¡¯¡@¡¯¡@¡¯¡@¡¯¡@¢x¡@¢x¡@¢x¡@¢x¡@¡¯¡@¡¯¡@¡¯¡@¡¯¡@¢x
        % ¢|¡C¢r¡C¢r¡C¢r¡C¢r¡@¢r¡@¢r¡@¢r¡@¢r¡C¢r¡C¢r¡C¢r¡C¢}
        %  ¡ö¡@ ¢Þ¢Û¢Ú ¡@¡÷¡@¡@¡@¡@¡@¡@¡@¡ö¡@ ¢Þ¢Û¢Ú ¡@¡÷
            %% Construct the optimal grids for PML layers
        % The grid of ele_left is determined by the bisection iteration.
        [ pml.optimal_factor, ele_left_grid, ele_left_mesh_length, ele_left_thickness, ele_left_thickness_ratio ] = ...
            FAME_Boundary_Perfect_Matching_Layer_Optimal_Mesh( grid_num, pml_grid_num, layer_thickness, mesh_len, pml_opt );
        % The grid of mag_left is the mid-point of ele_left.
        mag_left_grid = ( ele_left_grid(1:end-1) + ele_left_grid(2:end) ) / 2;
        % The grid of mag_right is symmetric to ele_left.
        mag_right_grid = edge_len - mesh_len/2 - flipud(ele_left_grid);
        % The grid of ele_right is the mid-point of mag_right.
        ele_right_grid = ( mag_right_grid(1:end-1) + mag_right_grid(2:end) ) / 2;
%         figure(1); 
%         hold on       
%         plot( mag_left_grid, ones(length( mag_left_grid),1), 'r*')
%         plot(mag_right_grid, ones(length(mag_right_grid),1), 'g*')
%         plot(ele_right_grid, ones(length(ele_right_grid),1), 'b*')
%         c = find( (mesh_len*(1:1000)>edge_len),1, 'first');
%         plot(mesh_len*(0:c-1), ones(c,1), 'k*')
        %% Compute thickness and ratio of PML ( ratio will be used in PML scale_factor )
        mag_left_thickness       = layer_thickness - mag_left_grid;
        mag_left_thickness_ratio = mag_left_thickness/layer_thickness;
        
        mag_right_thickness       = mag_right_grid(2:end) - mag_right_grid(1);
        mag_right_thickness_ratio = mag_right_thickness/layer_thickness;
        
        ele_right_thickness       = ele_right_grid - mag_right_grid(1);
        ele_right_thickness_ratio = ele_right_thickness/layer_thickness;
        
        %% Compute mesh_length 
        % The mesh_length of ele_left has been determined by the bisection iteration.
        mag_left_mesh_length  = [mag_left_grid(2:end); ele_left_grid(end)] - mag_left_grid;
        ele_right_mesh_length = flipud(mag_left_mesh_length);
        mag_right_mesh_length = flipud(ele_left_mesh_length);
        
        %% Construct the discrete nodes on electric and magnetic fields for perfect matching layer function s(x)
        % max_loss_parameter = -log(R)*(m+1)/(2*d_i*£s*£`_0*£b)
        max_loss_parameter = 8*( pml_opt.power+1 ) / ( pml_opt.frequency*layer_thickness );
        ele_left_scale_factor  = 1./( 1 - 1i*max_loss_parameter*ele_left_thickness_ratio .^pml_opt.power );
        mag_left_scale_factor  = 1./( 1 - 1i*max_loss_parameter*mag_left_thickness_ratio .^pml_opt.power );
        ele_right_scale_factor = 1./( 1 - 1i*max_loss_parameter*ele_right_thickness_ratio.^pml_opt.power );
        mag_right_scale_factor = 1./( 1 - 1i*max_loss_parameter*mag_right_thickness_ratio.^pml_opt.power );
        
        %% Constructing Boundary matrices ele_i and mag_i
        scale.ele_left_matrix  = [                  ele_left_scale_factor ; inner_one ; ele_right_scale_factor                  ];
        scale.mag_left_matrix  = [                  mag_left_scale_factor ; inner_one ; mag_right_scale_factor                  ];
        
        numesh.ele_left_matrix = [ mesh_len./ele_left_mesh_length ; inner_one ; mesh_len./ele_right_mesh_length ];
        numesh.mag_left_matrix = [ mesh_len./mag_left_mesh_length ; inner_one ; mesh_len./mag_right_mesh_length ];
        %%
        pml.ele_left_matrix      =    scale.ele_left_matrix.*numesh.ele_left_matrix;
        pml.mag_left_matrix      =    scale.mag_left_matrix.*numesh.mag_left_matrix;
        pml.ele_left_matrix_inv  = 1./pml.ele_left_matrix;
        pml.mag_left_matrix_inv  = 1./pml.mag_left_matrix;
%         pml.ele_left_lowrank     = 0.5*sparse( grid_num*ones(grid_num,1), 1:grid_num, pml.ele_left_matrix,...
%                                                            grid_num, grid_num );
%         pml.ele_left_lowrank_inv =     sparse( grid_num*ones(grid_num,1), 1:grid_num, pml.ele_left_matrix_inv(end)*ones(grid_num,1),...
%                                                            grid_num, grid_num );
% 
%         pml.mag_right_matrix     = speye(grid_num) - 0.5*sparse( 1:grid_num, grid_num*ones(grid_num,1), ones(grid_num,1) );
%         pml.mag_right_matrix_inv = speye(grid_num) +     sparse( 1:grid_num, grid_num*ones(grid_num,1), ones(grid_num,1) );

%         % Constructing the discrete nodes on electric field for perfect matching layer function s(x) 
%         ele_rho = linspace(layer_thickness, mesh_len, pml_grid_number).';
%         % Constructing the discrete nodes on magnetic field for perfect matching layer function s(x) 
%         mag_rho = ele_rho - (mesh_len / 2);
%         % Constructing the perfect matching layer function s(x) on electric field
%         ele_s = log(1/pml_opt.ratio)*((pml_opt.power+1)/(2*pml_opt.frequency*layer_thickness))*...
%                                (ele_rho./layer_thickness).^pml_opt.power;
%         % Constructing the perfect matching layer function s(x) on magnetic field
%         mag_s = log(1/pml_opt.ratio)*((pml_opt.power+1)/(2*pml_opt.frequency*layer_thickness))*...
%                                (mag_rho./layer_thickness).^pml_opt.power;
%         % Constructing the perfect matching layer function 1/(1-i*s(x)) on electric field       
%         ele_layer = 1./(1-1i*ele_s);
%         % Constructing the perfect matching layer function 1/(1-i*s(x)) on magnetic field        
%         mag_layer = 1./(1-1i*mag_s);
% 
%         %% Constructing Boundary matrices ele_i and mag_i
%         pml_ele_left_matrix = [ ele_layer ; inner_one ; flipud( ele_layer ) ];
%         pml_mag_left_matrix = [ mag_layer ; inner_one ; flipud( mag_layer ) ];
    else
        warning('grid_number must larger than 2*pml_grid_number!')
    end

end