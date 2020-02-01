function Call_Back_lattice_constant(hobj,event,h)
    
    [Popt] = FAME_GUI_Info_Get(h);

    hobj_idx   = get(hobj,'value');
    hobj_value = str2num(get(hobj,'string'));
    
    lattice_constant = Lattice_Constants_Pretreatment(Popt.lattice.lattice_type,Popt.lattice.lattice_constant,hobj_idx,hobj_value);
    
    grid_nums = Popt.mesh.grid_num;

    Pmaterial.data_name       = Popt.material.data_name;
    Pmaterial.sphere_radius   = Popt.material.sphere_radius;
    Pmaterial.cylinder_radius = Popt.material.cylinder_radius;

%     [ lattice_vec_a, edge_lens, mesh_lens, Omega, lattice_type, lattice_constant, Pmaterial ] = ...
%         FAME_Parameter_Generator( grid_nums, Popt.lattice.lattice_type, lattice_constant, Pmaterial );
    [ Par_mesh, Par_lattice ] = FAME_Parameter_Generator( grid_nums, Popt.lattice.lattice_type, lattice_constant, Pmaterial);
    lattice_vec_a = Par_lattice.Omega\Par_lattice.lattice_vec_a;
    lattice_constant = Par_lattice.lattice_constant;
    edge_lens = Par_mesh.edge_len;
    mesh_lens = Par_mesh.mesh_len;
    lattice_type = Par_lattice.lattice_type;
    
    cla(h.ax_show3D_prim)
    cla(h.ax_show3D_unit)
    FAME_Main_User_Defined_gui(Popt.material.data_name,Popt.material.sphere_radius,Popt.material.cylinder_radius, ...
        lattice_type, lattice_constant, lattice_vec_a, h.ax_show3D_prim, h.ax_show3D_unit);

    lattice_constant.a     = num2str(lattice_constant.a);
    lattice_constant.b     = num2str(lattice_constant.b);
    lattice_constant.c     = num2str(lattice_constant.c);
    lattice_constant.alpha = num2str(lattice_constant.alpha);
    lattice_constant.beta  = num2str(lattice_constant.beta);
    lattice_constant.gamma = num2str(lattice_constant.gamma);
    lattice_vec = {num2str(lattice_vec_a(:,1)), num2str(lattice_vec_a(:,2)), num2str(lattice_vec_a(:,3))};
    edge_len    = {num2str(edge_lens(1)), num2str(edge_lens(2)), num2str(edge_lens(3))};
    mesh_len    = {num2str(mesh_lens(1)), num2str(mesh_lens(2)), num2str(mesh_lens(3))};
    
    FAME_GUI_Info_Set('lattice_vector', lattice_vec, h);
    FAME_GUI_Info_Set('edge_length', edge_len, h);
    FAME_GUI_Info_Set('mesh_length', mesh_len, h);
    FAME_GUI_Info_Set('lattice_type', lattice_type, h);
    FAME_GUI_Info_Set('lattice_constant', lattice_constant, h);
    
    set( h.edit_domain.lattice_constant(1), 'string', lattice_constant.a);
    set( h.edit_domain.lattice_constant(2), 'string', lattice_constant.b);
    set( h.edit_domain.lattice_constant(3), 'string', lattice_constant.c);
    set( h.edit_domain.lattice_constant(4), 'string', lattice_constant.alpha);
    set( h.edit_domain.lattice_constant(5), 'string', lattice_constant.beta);
    set( h.edit_domain.lattice_constant(6), 'string', lattice_constant.gamma);
end

function lattice_constant = Lattice_Constants_Pretreatment(lattice_type,lattice_constant,hobj_idx,hobj_value)
% This routine reture pretreatment lattice constant which satisfying
%               a >= b >= c
    lattice_constant.old = lattice_constant;
    switch lattice_type
        case {'simple_cubic','face_centered_cubic','body_centered_cubic'}
            if ismember(hobj_idx,[1 2 3])
                lattice_constant.a     = hobj_value;
                lattice_constant.b     = hobj_value;
                lattice_constant.c     = hobj_value;
            end
            lattice_constant.alpha = pi/2;
            lattice_constant.beta  = pi/2;
            lattice_constant.gamma = pi/2;
        case {'hexagonal' ,'rhombohedral'}
            if ismember(hobj_idx,[1 2])
                lattice_constant.a     = hobj_value;
                lattice_constant.b     = hobj_value;
            elseif hobj_idx == 3
                lattice_constant.c     = hobj_value;
            end
            lattice_constant.alpha = pi/2;
            lattice_constant.beta  = pi/2;
            lattice_constant.gamma = 2*pi/3;
        case {'primitive_tetragonal','body_centered_tetragonal'}
            if ismember(hobj_idx,[1 2])
                lattice_constant.a     = hobj_value;
                lattice_constant.b     = hobj_value;
            elseif hobj_idx == 3
                lattice_constant.c     = hobj_value;
            end
            lattice_constant.alpha = pi/2;
            lattice_constant.beta  = pi/2;
            lattice_constant.gamma = pi/2;
        case {'primitive_orthorhombic'}
            if hobj_idx == 1
                lattice_constant.a     = hobj_value;
            elseif hobj_idx == 2    
                lattice_constant.b     = hobj_value;
            elseif hobj_idx == 3
                lattice_constant.c     = hobj_value;
            end
            lattice_constant.alpha = pi/2;
            lattice_constant.beta  = pi/2;
            lattice_constant.gamma = pi/2;
        case {'a_base_centered_orthorhombic','c_base_centered_orthorhombic',...
               'face_centered_orthorhombic', 'body_centered_orthorhombic'}
            if hobj_idx == 1
                lattice_constant.a     = hobj_value;
            elseif hobj_idx == 2    
                lattice_constant.b     = hobj_value;
            elseif hobj_idx == 3
                lattice_constant.c     = hobj_value;
            end
            lattice_constant.alpha = pi/2;
            lattice_constant.beta  = pi/2;
            lattice_constant.gamma = pi/2;    
        case {'primitive_monoclinic','a_base_centered_monoclinic'}
            if hobj_idx == 1
                lattice_constant.a     = hobj_value;
            elseif hobj_idx == 2    
                lattice_constant.b     = hobj_value;
            elseif hobj_idx == 3
                lattice_constant.c     = hobj_value;
            elseif hobj_idx == 4
                lattice_constant.alpha = hobj_value;
                lattice_constant.beta  = pi/2;
                lattice_constant.gamma = pi/2;
            elseif hobj_idx == 5
                lattice_constant.alpha = pi/2;
                lattice_constant.beta  = hobj_value;
                lattice_constant.gamma = pi/2;
            elseif hobj_idx == 6
                lattice_constant.alpha = pi/2;    
                lattice_constant.beta  = pi/2;
                lattice_constant.gamma = hobj_value;
            end
        case 'triclinic'    
            if     hobj_idx == 1
                lattice_constant.a     = hobj_value;
            elseif hobj_idx == 2    
                lattice_constant.b     = hobj_value;
            elseif hobj_idx == 3
                lattice_constant.c     = hobj_value;
            elseif hobj_idx == 4
                lattice_constant.alpha = hobj_value;
            elseif hobj_idx == 5
                lattice_constant.beta  = hobj_value;
            elseif hobj_idx == 6
                lattice_constant.gamma = hobj_value;
            end
    end
%     lattice_constant.b = lattice_constant.b/lattice_constant.a;
%     lattice_constant.c = lattice_constant.c/lattice_constant.a;
%     lattice_constant.a = 1;
end
