function Par_recip_lattice = FAME_Parameter_Brillouin_Zone_Path( varargin )
    flag_usr = 0;
    switch nargin 
        case {0,1}
            error('To less input argument!')
        case 2
            wave_vec_array_usr = varargin{1};
            Par_recip_lattice  = varargin{2};
            flag_usr = 1;
            path_string = [];
            path_string_new = [];
            vertex = [];
            part_num = [];
        case 3
            part_num          = varargin{1};
            Par_lattice       = varargin{2};
            Par_recip_lattice = varargin{3};
            path_string       = default_path(Par_lattice.lattice_type, Par_lattice.lattice_constant);
        case 4
            part_num          = varargin{1};
            Par_lattice       = varargin{2};
            Par_recip_lattice = varargin{3};
            path_string       = varargin{4};
    end
    
    if flag_usr == 1
        wave_vec_array_usr = reshape(wave_vec_array_usr,3,size(wave_vec_array_usr,1)*size(wave_vec_array_usr,2)/3);
        wave_vec_array     = Par_recip_lattice.reciprocal_lattice_vector_b*wave_vec_array_usr;
    else
        [ vertex_orig ] = FAME_Parameter_Brillouin_Zone_Point(Par_lattice.lattice_type, Par_lattice.lattice_constant_orig, Par_recip_lattice.reciprocal_lattice_vector_b_orig);
        vertex_str = fieldnames(vertex_orig);
        vertex = vertex_orig;
        for i = 1:length(vertex_str)
            eval( ['vertex.',vertex_str{i},' = Par_lattice.Omega*vertex_orig.',vertex_str{i},';'] );
        end
        
        n           = length(path_string);
        wave_vec_array  = [];
        count_part  = 0;
        path_string = [ path_string, '-'];
        n_new       = 1;
        flag        = 0;
        if n == 1 
            eval( ['wave_vec_array = [ wave_vec_array, vertex.',path_string(1),'];']);
        end
        for i = 1:n-1
            if strcmp(path_string(i),'|') == 0 && strcmp(path_string(i+1),'|') == 0
                subpath_start_string = [ 'vertex.', path_string(i)   ];
                subpath_end_string   = [ 'vertex.', path_string(i+1) ];
                string_x = [ 'subpath(1,:) = linspace(',subpath_start_string,'(1),',subpath_end_string,'(1),', num2str(part_num),');'];
                string_y = [ 'subpath(2,:) = linspace(',subpath_start_string,'(2),',subpath_end_string,'(2),', num2str(part_num),');'];
                string_z = [ 'subpath(3,:) = linspace(',subpath_start_string,'(3),',subpath_end_string,'(3),', num2str(part_num),');'];

                eval(string_x);
                eval(string_y);
                eval(string_z);
                if i == n-1 || strcmp(path_string(i+2),'|') == 1
                    wave_vec_array = [wave_vec_array,subpath(:,1:end)];
                else
                    wave_vec_array = [wave_vec_array,subpath(:,1:end-1)];
                end
            end
            if strcmp(path_string(i),'|') == 1
                count_part = count_part + 1;
                path_string_new{n_new} = path_string(i-1:i+1);
                n_new = n_new + 1;
                flag  = 1;
            elseif strcmp(path_string(i+1),'|') == 1
            elseif flag == 1
                flag = 0;
            else
                path_string_new{n_new} = path_string(i);
                n_new = n_new + 1;
            end
        end
        path_string_new{n_new} = path_string(i+1);
    end
    
    wave_vec_num           = size(wave_vec_array,2);

    Par_recip_lattice.wave_vec_array = wave_vec_array;
    Par_recip_lattice.wave_vec_num = wave_vec_num; 
    Par_recip_lattice.path_string = path_string;
    Par_recip_lattice.path_string_new = path_string_new;
    Par_recip_lattice.vertex = vertex; 
    Par_recip_lattice.part_num = part_num;
end

    
function path_string = default_path(lattice_type, lattice_constant)
    switch lattice_type
        case 'simple_cubic'   
            path_string = 'GXMGR|MR';
        case 'face_centered_cubic'
            path_string = 'GXWKGLUWLK|UX';
        case 'body_centered_cubic'
            path_string = 'GHNGPH|PN';
        case 'hexagonal'    
            path_string = 'GMKGALHA|LM|KH';
        case 'rhombohedral'
            theta = 2*asin( 3*lattice_constant.a/(2*sqrt( lattice_constant.c^2 + 3*lattice_constant.a^2 ) ));
            if theta <= pi/2
            path_string = 'GLC|BZGX|QFRZ|LP';
            elseif theta > pi/2
            path_string = 'GPZQGFRSLZ';
            end
        case 'primitive_tetragonal'
            path_string = 'GXMGZRAZ|XR|MA';
        case 'body_centered_tetragonal'
            if lattice_constant.c < lattice_constant.a
            path_string = 'GXMGZPNBM|XP';
            elseif lattice_constant.c > lattice_constant.a
            path_string = 'GXYsGZtNPAZ|XP';
            end
        case 'primitive_orthorhombic'
%             path_string = 'GXSYGZURTZ|YT|UX|SR';
%             path_string = 'RYSRTYGXUGZUR';
            path_string = 'RXSRSXGZUGZUR';
    %     case 'a_base_centered_orthorhombic'
    %         path_string = '';
        case 'c_base_centered_orthorhombic'
            path_string = 'GXSRAZGYCBTY|ZT';
        case 'face_centered_orthorhombic'
            if 1/lattice_constant.a^2 > 1/lattice_constant.b^2 + 1/lattice_constant.c^2
            path_string = 'GYTZGXBY|TJ|XAZ|LG';
            elseif 1/lattice_constant.a^2 < 1/lattice_constant.b^2 + 1/lattice_constant.c^2
            path_string = 'GYCDXGZFHC|EZ|XI|HY|LG';
            elseif 1/lattice_constant.a^2 == 1/lattice_constant.b^2 + 1/lattice_constant.c^2
            path_string = 'GYTZGXBY|XAZ|LG';
            end
        case 'body_centered_orthorhombic'
            path_string = 'GXLTWRAZGYSW|IY|BZ';
        case 'primitive_monoclinic'
            path_string = 'GYHCENAXI|MDZ|YD';
        case 'base_centered_monoclinic'
            a = lattice_vec_a;
            b = inv(a)';
            temp = lattice_constant.b  *cos(lattice_constant.alpha)  /lattice_constant.c   + ...
                   lattice_constant.b^2*sin(lattice_constant.alpha)^2/lattice_constant.a^2;     

            k_gamma = acos( (b(:,1).'*b(:,2))/(norm(b(:,1))*norm(b(:,2))) );
            if k_gamma > pi/2
                path_string = 'GYFLI|RZJ|YT|XGN|MG';
            elseif k_gamma == pi/2
                path_string = 'GYFLI|RZJ|NGM';
            elseif k_gamma < pi/2 && temp < 1
                path_string = 'GYFHZIJ|PVXGN|MG';
            elseif k_gamma < pi/2 && temp == 1
                path_string = 'GYFHZI|PVXGN|MG';
            elseif k_gamma < pi/2 && temp > 1
                path_string = 'GYFLI|RZHJ|PVXGN|MG';
            end
        case 'triclinic'    
            path_string = 'XGY|LGZ|NGM|RG';
    end
end