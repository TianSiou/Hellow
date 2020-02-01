function varargout = FAME_Main_Code_nonorth( varargin )
%FAME_Main_Code Compute the source-free Maxwell's equations of a photonic
%crystal. 

% warning off

[ Par_mesh, Par_lattice, Par_recip_lattice, Par_material, Par_eig, wave_vec_array, FAME_option ] = input_check(nargin,nargout,varargin);

fprintf('--------------------------------Start FAME program----------------------------------\n')
%% (Main Step 1) - Construct the matrices B's
fprintf('Construct material dependent matrix B......')
switch FAME_option.discrete_method
    case 'Yee_scheme'
        switch Par_material.material_type
            case 'isotropic'
            B = FAME_Matrix_B_Isotropic_nonorth( Par_mesh, Par_material);
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
    edge_len = [1,1,1];
    mesh_len = [1/Par_mesh.grid_num(1), 1/Par_mesh.grid_num(2), 1/Par_mesh.grid_num(3)];
    SC_lattice_vec_a = [1,0,0;0,1,0;0,0,1];
    SC_lattice_constant.a = 1;
    SC_lattice_constant.b = 1;
    SC_lattice_constant.c = 1;
    SC_lattice_constant.alpha = pi/2;
    SC_lattice_constant.beta  = pi/2;
    SC_lattice_constant.gamma = pi/2;
    [ C , Cs, C_1, C_2, C_3, C_x, C_y, C_z ] = ...
        FAME_Matrix_Curl_nonorth( wave_vec, Par_mesh.grid_num, edge_len, mesh_len, {'quasi_periodic','quasi_periodic','quasi_periodic'}, 'simple_cubic', SC_lattice_vec_a, SC_lattice_constant );
    [ Lambdas ] = ...
        FAME_Matrix_Lambdas( wave_vec, Par_mesh.grid_num, mesh_len, 'simple_cubic', SC_lattice_constant, SC_lattice_vec_a );
    
    FAME_Check_Eigendecomp(Par_mesh.grid_num, 'simple_cubic', C, Cs, C_1, C_2, C_3, Lambdas);
    
    lattice_vec_b = Par_recip_lattice.reciprocal_lattice_vector_b;
    Ce_1 =  ( lattice_vec_b(1,1)*C_1  + lattice_vec_b(1,2)*C_2  + lattice_vec_b(1,3)*C_3  ) / (2*pi);
    Ce_2 =  ( lattice_vec_b(2,1)*C_1  + lattice_vec_b(2,2)*C_2  + lattice_vec_b(2,3)*C_3  ) / (2*pi);
    Ce_3 =  ( lattice_vec_b(3,1)*C_1  + lattice_vec_b(3,2)*C_2  + lattice_vec_b(3,3)*C_3  ) / (2*pi);
    Ch_1 = -( lattice_vec_b(1,1)*C_1' + lattice_vec_b(1,2)*C_2' + lattice_vec_b(1,3)*C_3' ) / (2*pi);
    Ch_2 = -( lattice_vec_b(2,1)*C_1' + lattice_vec_b(2,2)*C_2' + lattice_vec_b(2,3)*C_3' ) / (2*pi);
    Ch_3 = -( lattice_vec_b(3,1)*C_1' + lattice_vec_b(3,2)*C_2' + lattice_vec_b(3,3)*C_3' ) / (2*pi);
    
    O = sparse(size(Ce_1,1),size(Ce_1,2));
    Ce = [     O, -Ce_3,  Ce_2;
            Ce_3,     O, -Ce_1;
           -Ce_2,  Ce_1,     O];
    Ch = [     O, -Ch_3,  Ch_2;
            Ch_3,     O, -Ch_1;
           -Ch_2,  Ch_1,     O];
    EW = eig(full(Ch*Ce),diag(B.B_eps));
    EW_tmp = EW(EW>1e-6); EW_tmp = sort(EW_tmp,'ascend'); EW_tmp = EW_tmp(1:6);
    Freq_array(:,i) = sqrt(EW_tmp);
    if norm(wave_vec) < (1e-6)
        Freq_array(:,i) = [0;0;Freq_array(1:end-2,i)];
    end
    
    figure(3);
    plot(i*ones(6,1),EW_tmp,'b.'); hold on
    
    clear('Lambdas','wave_vec','C','Cs');
    fprintf('================================== %3.1f%% complete =====================================\n',100*(i/size(wave_vec_array,2)))
end

varargout = output_check(nargout, Freq_array);

end

function [ Par_mesh, Par_lattice, Par_recip_lattice, Par_material, Par_eig, Par_recip_lattice_wave_vec_array, FAME_option ] = input_check(Nargin,Nargout,Varargin)
    
    if Nargout > 4
        error('Too many output argument.');
    end
    switch Nargin
        case {0,1,2,3,4,5} 
           error('Too less input argument.')
        case 6
            Par_mesh                         = Varargin{1};
            Par_lattice                      = Varargin{2};
            Par_recip_lattice                = Varargin{3};
            Par_material                     = Varargin{4};
            Par_eig                          = Varargin{5};
            Par_recip_lattice_wave_vec_array = Varargin{6};
            FAME_option.discrete_method      = 'Yee_scheme';
        case 7
            Par_mesh                         = Varargin{1};
            Par_lattice                      = Varargin{2};
            Par_recip_lattice                = Varargin{3};
            Par_material                     = Varargin{4};
            Par_eig                          = Varargin{5};
            Par_recip_lattice_wave_vec_array = Varargin{6};
            FAME_option                      = Varargin{7};
        otherwise
            error('Too many input argument.');
    end 
end

function Varargout = output_check(Nargout, Freq_array)
    switch Nargout
        case {0,1}
            Varargout{1} = Freq_array;
        otherwise
            error('Too many output argument.');
    end
end

function [ Freq_array, Ele_field_mtx, Mag_field_mtx, cpu_time, LS_iter, LS_cpu_time, Err ] = FAME_Fast_Algorithms_Driver_yeescheme( Par_material, Par_lattice, grid_num, wave_vec, B, Lambdas, eigen_wanted, C, Cs )
    switch Par_material.material_type
        case 'isotropic'
        Mag_field_mtx = [];
        switch Par_lattice.lattice_type
            case {'simple_cubic', 'primitive_orthorhombic', 'primitive_tetragonal'}
            [ Freq_array, Ele_field_mtx, cpu_time, LS_iter, LS_cpu_time ] = ...
                FAME_Fast_Algorithms_Simple_Isotropic( grid_num, wave_vec, B, Lambdas, eigen_wanted);
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