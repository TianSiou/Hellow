function varargout = FAME_Main_Code( varargin )
%FAME_Main_Code Compute the source-free Maxwell's equations of a photonic
%crystal. 

% warning off
[ Par_mesh, Par_lattice, Par_material, Par_eig, wave_vec_array, FAME_option ] = input_check(nargin,nargout,varargin);

fprintf('--------------------------------Start FAME program----------------------------------\n')
%% (Main Step 1) - Construct the matrices B's
fprintf('Construct material dependent matrix B......')
switch FAME_option.discrete_method
    case 'Yee_scheme'
        switch Par_material.material_type
            case 'isotropic'
            B = FAME_Matrix_B_Isotropic( Par_mesh, Par_material);
            case 'anisotropic'
            B = FAME_Matrix_B_Anisotropic( Par_mesh, Par_material );
            case 'biisotropic'
            B = FAME_Matrix_B_Biisotropic( Par_mesh, Par_material );    
        end
    case 'fem'
        switch Par_material.material_type
            case 'isotropic'
                % Not valid now
            case 'anisotropic'
            B = FAME_Matrix_B_Anisotropic_fem_1( Par_mesh, Par_lattice, Par_material);
            case 'biisotropic'
            B = FAME_Matrix_B_Biisotropic_fem_1( Par_mesh, Par_lattice, Par_material );    
        end
end
fprintf('Done\n')
% Display compuational information
FAME_Tools_Display_Info( Par_mesh, Par_lattice, Par_material, wave_vec_array )
fprintf('===========================Start computing band structure===========================\n')
for i = 1:size(wave_vec_array,2)
    % Choose wave vector in Brillouin zone path
    wave_vec = wave_vec_array(:,i);
    %% (Main Step 2) - Construct Lambdas and discrete curl-operators corresponding to the choosen wave vector
    fprintf('Construct Lambdas...')
    [ C , Cs, C_1, C_2, C_3, C_x, C_y, C_z ] = ...
        FAME_Matrix_Curl( wave_vec, Par_mesh.grid_num, Par_mesh.edge_len,Par_mesh.mesh_len, {'quasi_periodic','quasi_periodic','quasi_periodic'}, Par_lattice.lattice_type, Par_lattice.lattice_vec_a, Par_lattice.lattice_constant );
%     [ C_sc , Cs_sc, C_1_sc, C_2_sc, C_3_sc, C_x_sc, C_y_sc, C_z_sc ] = ...
%         FAME_Matrix_Curl( wave_vec, [16, 16, 16], [1, 1, 1], Par_mesh.mesh_len, {'quasi_periodic','quasi_periodic','quasi_periodic'}, Par_lattice.lattice_type, eye(3), Par_lattice.lattice_constant );
    
    [ Lambdas ] = ...
        FAME_Matrix_Lambdas( wave_vec, Par_mesh.grid_num, Par_mesh.mesh_len, Par_lattice.lattice_type, Par_lattice.lattice_constant, Par_lattice.lattice_vec_a );
    
    FAME_Check_Eigendecomp(Par_mesh.grid_num, Par_lattice.lattice_type, C, Cs, C_1, C_2, C_3, Lambdas);
    
%     shift =  pi^2;
%     N = size(C,1);
%    
%     [ev, ew, flag] = eigs( @(x) shift_inv_BCsC(x, B.invB_eps, Cs, C, shift), N, 1, 'lm');
%     ew = diag(ew);
%     ew = sqrt(ew)/2/pi;
    fprintf('Done\n')
    if strcmp(FAME_option.discrete_method,'fem') == 1
        switch Par_material.material_type
            case 'anisotropic'
                B = FAME_Matrix_B_Anisotropic_fem_2( B, C_1, C_2, C_3, Par_mesh.mesh_len(1), Par_mesh.mesh_len(2), Par_mesh.mesh_len(3));
            case 'biisotropic'
                B = FAME_Matrix_B_Biisotropic_fem_2( B, C_1, C_2, C_3, Par_mesh.mesh_len(1), Par_mesh.mesh_len(2), Par_mesh.mesh_len(3));
        end
    end
    
   Info.EP_cpu_time = [];
    %% (Main Step 3) - Use fast algorithms to compute the smallest frequency for given wave vector, material type, and lattice type
    fprintf('Start FAME(%s, wave vector = [%.2e,%.2e,%.2e])...', Par_material.material_type, wave_vec)
    switch FAME_option.discrete_method
        case 'Yee_scheme'
            [ Freq_array(:,i), Ele_field_mtx{i}, Mag_field_mtx{i}, Info.cpu_time(i), Info.LS_iter{i}, Info.LS_cpu_time{i}, Info.Err(:,i) ,Info.EP_cpu_time{i}, Info.Sira_iter{i}] = ...
                FAME_Fast_Algorithms_Driver_yeescheme( Par_material, Par_lattice, Par_mesh.grid_num, wave_vec, B, Lambdas, Par_eig.eigen_wanted, C, Cs );
%             [ Freq_array(:,i), Ele_field_mtx{i}, Mag_field_mtx{i}, Info.cpu_time(i), Info.LS_iter{i}, Info.LS_cpu_time{i}, Info.Err(:,i) ] = ...
%                 FAME_Fast_Algorithms_Driver_yeescheme( Par_material, Par_lattice, Par_mesh.grid_num, wave_vec, B, Lambdas, Par_eig.eigen_wanted, C, Cs , C_sc, Cs_sc);
        case 'fem'
            [ Freq_array(:,i), Ele_field_mtx{i}, Mag_field_mtx{i}, Info.cpu_time(i), Info.LS_iter{i}, Info.LS_cpu_time{i}, Info.Err(:,i) ] = ...
                FAME_Fast_Algorithms_Driver_fem( Par_material, Par_lattice, Par_mesh.grid_num, wave_vec, B, Lambdas, Par_eig.eigen_wanted, C, Cs, C_x, C_y, C_z );
    end
    
    fprintf('Done! Use time %.2f(sec.)\n',Info.cpu_time(i))
    
    clear('Lambdas','wave_vec','C','Cs');
    fprintf('================================== %3.1f%% complete =====================================\n',100*(i/size(wave_vec_array,2)))
end

varargout = output_check(nargout, Freq_array, Ele_field_mtx, Mag_field_mtx,Info);

end

function vec_y = shift_inv_BCsC(vec_x, invB_eps, Cs, C, shift)
 [ vec_y, flag, res, inner_iter] = bicgstab(  @(x) invB_eps .*( Cs * (C * x) ) - shift * x, vec_x, 1e-8, 10000);
inner_iter
end


function [ Par_mesh, Par_lattice, Par_material, Par_eig, Par_recip_lattice_wave_vec_array, FAME_option ] = input_check(Nargin,Nargout,Varargin)
    
    if Nargout > 4
        error('Too many output argument.');
    end
    switch Nargin
        case {0,1,2,3,4} 
           error('Too less input argument.')
        case 5
            Par_mesh                         = Varargin{1};
            Par_lattice                      = Varargin{2};
            Par_material                     = Varargin{3};
            Par_eig                          = Varargin{4};
            Par_recip_lattice_wave_vec_array = Varargin{5};
            FAME_option.discrete_method      = 'Yee_scheme';
        case 6
            Par_mesh                         = Varargin{1};
            Par_lattice                      = Varargin{2};
            Par_material                     = Varargin{3};
            Par_eig                          = Varargin{4};
            Par_recip_lattice_wave_vec_array = Varargin{5};
            FAME_option                      = Varargin{6};
        otherwise
            error('Too many input argument.');
    end 
end

function Varargout = output_check(Nargout, Freq_array, Ele_field_mtx, Mag_field_mtx, Info)
    switch Nargout
        case {0,1,2}
            error('Too less output argument.')
        case 3
            Varargout{1} = Freq_array;
            Varargout{2} = Ele_field_mtx;
            Varargout{3} = Mag_field_mtx;
        case 4
            Varargout{1} = Freq_array;
            Varargout{2} = Ele_field_mtx;
            Varargout{3} = Mag_field_mtx;    
            Varargout{4} = Info;    
        otherwise
            error('Too many output argument.');
    end
end

function [ Freq_array, Ele_field_mtx, Mag_field_mtx, cpu_time, LS_iter, LS_cpu_time, Err, EP_cpu_time, Sira_iter] = FAME_Fast_Algorithms_Driver_yeescheme( Par_material, Par_lattice, grid_num, wave_vec, B, Lambdas, eigen_wanted, C, Cs)
    switch Par_material.material_type
        case 'isotropic'
        EP_cpu_time = 0;
        Mag_field_mtx = [];
        switch Par_lattice.lattice_type
            case {'simple_cubic', 'primitive_orthorhombic', 'primitive_tetragonal'}
            [ Freq_array, Ele_field_mtx, cpu_time, LS_iter, LS_cpu_time , EP_cpu_time, Sira_iter] = ...
                FAME_Fast_Algorithms_Simple_Isotropic_test( grid_num, wave_vec, B, Lambdas, eigen_wanted);         
%             [ Freq_array, Ele_field_mtx, cpu_time, LS_iter, LS_cpu_time ] = ...
%                 FAME_Fast_Algorithms_Simple_Isotropic( grid_num, wave_vec, B, Lambdas, eigen_wanted);
%             [ Freq_array, Ele_field_mtx, cpu_time, LS_iter, LS_cpu_time ] = ...
%                 FAME_Fast_Algorithms_Simple_Isotropic_gpu( grid_num, wave_vec, B, Lambdas, eigen_wanted);
            otherwise
            [ Freq_array, Ele_field_mtx, cpu_time, LS_iter, LS_cpu_time] = ...
                FAME_Fast_Algorithms_General_Isotropic(grid_num, wave_vec, B, Lambdas, eigen_wanted);
        end
        Err = FAME_Error_Check_Isotropic( Freq_array, Ele_field_mtx,C, Cs, B );
        
        case 'anisotropic'
        Mag_field_mtx = [];
        switch Par_lattice.lattice_type
            case {'simple_cubic', 'primitive_orthorhombic', 'primitive_tetragonal'}
            [ Freq_array, Ele_field_mtx, cpu_time, LS_iter, LS_cpu_time ] = ...
                FAME_Fast_Algorithms_Simple_Anisotropic(grid_num, wave_vec, B, Lambdas, eigen_wanted);
            otherwise
            [ Freq_array, Ele_field_mtx, cpu_time, LS_iter, LS_cpu_time ] = ...
                FAME_Fast_Algorithms_General_Anisotropic(grid_num, wave_vec, B, Lambdas, eigen_wanted);
        end
        Err = FAME_Error_Check_Anisotropic( Freq_array, Ele_field_mtx,C, Cs, B );
        case 'biisotropic'
        switch Par_lattice.lattice_type
            case {'simple_cubic', 'primitive_orthorhombic', 'primitive_tetragonal'}
            [ Freq_array, Ele_field_mtx, Mag_field_mtx, cpu_time, LS_iter, LS_cpu_time ] = ...
                FAME_Fast_Algorithms_Simple_Biisotropic(grid_num, wave_vec, B, Lambdas, eigen_wanted);
            otherwise
            [ Freq_array, Ele_field_mtx, Mag_field_mtx, cpu_time, LS_iter, LS_cpu_time ] = ...
                FAME_Fast_Algorithms_General_Biisotropic(grid_num, wave_vec, B, Lambdas, eigen_wanted);
        end
        Err = FAME_Error_Check_Biisotropic( Freq_array, Ele_field_mtx, Mag_field_mtx,C, Cs, B );
    end 
end

function [ Freq_array, Ele_field_mtx, Mag_field_mtx, cpu_time, LS_iter, LS_cpu_time, Err ] = FAME_Fast_Algorithms_Driver_fem( Par_material, Par_lattice, grid_num, wave_vec, B, Lambdas, eigen_wanted, C, Cs, C_x, C_y, C_z )
    switch Par_material.material_type
        case 'isotropic'
        Mag_field_mtx = [];
        switch Par_lattice.lattice_type
            case {'simple_cubic', 'primitive_orthorhombic', 'primitive_tetragonal'}
            [ Freq_array, Ele_field_mtx, cpu_time, LS_iter, LS_cpu_time ] = ...
                FAME_Fast_Algorithms_Simple_Isotropic_fem( grid_num, wave_vec, B, Lambdas, eigen_wanted);
            otherwise
            [ Freq_array, Ele_field_mtx, cpu_time, LS_iter, LS_cpu_time] = ...
                FAME_Fast_Algorithms_General_Isotropic_fem(grid_num, wave_vec, B, Lambdas, eigen_wanted);
        end
        Err = FAME_Error_Check_Isotropic_fem( Freq_array, Ele_field_mtx,C, Cs, B );
        
        case 'anisotropic'
        Mag_field_mtx = [];
        switch Par_lattice.lattice_type
            case {'simple_cubic', 'primitive_orthorhombic', 'primitive_tetragonal'}
            [ Freq_array, Ele_field_mtx, cpu_time, LS_iter, LS_cpu_time ] = ...
                FAME_Fast_Algorithms_Simple_Anisotropic_fem(grid_num, wave_vec, B, Lambdas, eigen_wanted);
            otherwise
            [ Freq_array, Ele_field_mtx, cpu_time, LS_iter, LS_cpu_time ] = ...
                FAME_Fast_Algorithms_General_Anisotropic_fem(grid_num, wave_vec, B, Lambdas, eigen_wanted);
        end
        Err = FAME_Error_Check_Anisotropic_fem( Freq_array, Ele_field_mtx,C, Cs, B );
        case 'biisotropic'
        switch Par_lattice.lattice_type
            case {'simple_cubic', 'primitive_orthorhombic', 'primitive_tetragonal'}
            [ Freq_array, Ele_field_mtx, Mag_field_mtx, cpu_time, LS_iter, LS_cpu_time ] = ...
                FAME_Fast_Algorithms_Simple_Biisotropic_fem(grid_num, wave_vec, B, Lambdas, eigen_wanted);
            otherwise
            [ Freq_array, Ele_field_mtx, Mag_field_mtx, cpu_time, LS_iter, LS_cpu_time ] = ...
                FAME_Fast_Algorithms_General_Biisotropic_fem(grid_num, wave_vec, B, Lambdas, eigen_wanted);
        end
        Err = FAME_Error_Check_Biisotropic_fem( Freq_array, Ele_field_mtx, Mag_field_mtx,C, Cs, C_x, C_y, C_z, B );
    end
    
end