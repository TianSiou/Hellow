function Call_Back_shape(hobj,event,h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Get informations
% Read choosing file name
popup_num                     = get(h.popup_shape,'value');
file_name_list                = get(h.popup_shape,'string');
file_name                     = file_name_list(popup_num,:);
file_name(isspace(file_name)) = [];

material = FAME_Material_Locate_Parameter_Get(file_name); % Get material properties

[Popt] = FAME_GUI_Info_Get(h); % Get grid numbers from information figure
% Generate lattice and mesh informations
grid_num                  = [Popt.mesh.grid_num(1), Popt.mesh.grid_num(2), Popt.mesh.grid_num(3)];
Pmaterial.data_name       = file_name;
Pmaterial.sphere_radius   = Popt.material.sphere_radius;
Pmaterial.cylinder_radius = Popt.material.cylinder_radius;

% Close some lattice constants corresponding to lattice type
switch material.lattice_type
    case {'simple_cubic','face_centered_cubic','body_centered_cubic'}
        Idx_close = 2:6;
        Idx_open  = 1;
    case {'hexagonal' ,'rhombohedral'}
        Idx_close = [2,4,5,6];
        Idx_open  = [1,3];
    case {'primitive_tetragonal','body_centered_tetragonal'}    
        Idx_close = [2,4,5,6];
        Idx_open  = [1,3];
    case {'primitive_orthorhombic','a_base_centered_orthorhombic','c_base_centered_orthorhombic','face_centered_orthorhombic', 'body_centered_orthorhombic'}
        Idx_close = [4,5,6];
        Idx_open  = [1,2,3];
    case {'primitive_monoclinic','base_centered_monoclinic'}
        Idx_close = [5,6];
        Idx_open  = [1,2,3,4];
    case 'triclinic' 
        Idx_close = [];
        Idx_open  = [1,2,3,4,5,6];
end
% if isempty(Idx_close) ~= 1
    set(h.edit_domain.lattice_constant(Idx_close), 'Enable', 'off')
% end
set(h.edit_domain.lattice_constant(Idx_open), 'Enable', 'on')


% [ lattice_vec_a, edge_len, mesh_len, Omega, lattice_type, lattice_constant, Pmaterial ] = ...
%     FAME_Parameter_Generator( grid_num, material.lattice_type, material.lattice_constant, Pmaterial );
[ Par_mesh, Par_lattice, ~, ~ ] = FAME_Parameter_Generator( grid_num, material.lattice_type, material.lattice_constant, []);
edge_len      = Par_mesh.edge_len;
mesh_len      = Par_mesh.mesh_len;
lattice_type  = Par_lattice.lattice_type;
lattice_constant = Par_lattice.lattice_constant;

lattice_vec_a = Par_lattice.Omega\Par_lattice.lattice_vec_a;

%% Plot material in primitive and unit cell
cla(h.ax_show3D_prim)
cla(h.ax_show3D_unit)

FAME_Main_User_Defined_gui(file_name,Popt.material.sphere_radius,Popt.material.cylinder_radius, ...
    lattice_type, lattice_constant, lattice_vec_a, h.ax_show3D_prim, h.ax_show3D_unit);

%% Update Information figure
lattice_constant.a     = num2str(lattice_constant.a);
lattice_constant.b     = num2str(lattice_constant.b);
lattice_constant.c     = num2str(lattice_constant.c);
lattice_constant.alpha = num2str(lattice_constant.alpha);
lattice_constant.beta  = num2str(lattice_constant.beta);
lattice_constant.gamma = num2str(lattice_constant.gamma);
lattice_constant.Permutation = num2str(lattice_constant.Permutation);

lattice_vec = {num2str(lattice_vec_a(:,1)), num2str(lattice_vec_a(:,2)), num2str(lattice_vec_a(:,3))};
edge_len    = {num2str(edge_len(1)), num2str(edge_len(2)), num2str(edge_len(3))};
mesh_len    = {num2str(mesh_len(1)), num2str(mesh_len(2)), num2str(mesh_len(3))};

FAME_GUI_Info_Set('lattice_vector', lattice_vec, h);
FAME_GUI_Info_Set('edge_length', edge_len, h);
FAME_GUI_Info_Set('mesh_length', mesh_len, h);
FAME_GUI_Info_Set('lattice_type', lattice_type, h);
FAME_GUI_Info_Set('lattice_constant', lattice_constant, h);
FAME_GUI_Info_Set('file_name', file_name, h);

%% Update main FAME figure
set( h.edit_domain.lattice_constant(1), 'string', lattice_constant.a );
set( h.edit_domain.lattice_constant(2), 'string', lattice_constant.b );
set( h.edit_domain.lattice_constant(3), 'string', lattice_constant.c );
set( h.edit_domain.lattice_constant(4), 'string', lattice_constant.alpha );
set( h.edit_domain.lattice_constant(5), 'string', lattice_constant.beta );
set( h.edit_domain.lattice_constant(6), 'string', lattice_constant.gamma );
end
