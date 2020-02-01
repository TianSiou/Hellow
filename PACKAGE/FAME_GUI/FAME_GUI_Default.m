%% Material default option
file_name_list = get(h.popup_shape,'string');
file_name      = 'No221_SC_NbO';

for i = 1:size(file_name_list,1)
    if strcmp(file_name_list(i,1:length(file_name)),file_name) == 1
        popup_num = i;
    end
end
set(h.popup_shape, 'value', popup_num);

Pmaterial.data_name       = file_name;
Pmaterial.sphere_radius   = 0.12;
Pmaterial.cylinder_radius = 0.1;

material = FAME_Material_Locate_Parameter_Get(file_name); % Get material properties

material_type = 'isotropic';

permittivity = '13';
permeability = '1';
reciprocity  = '0';
chirality    = '0';

part_num     = '3';

%% Mesh default option
grid_num                 = [6, 6, 6];

%% Lattice default option
% [ lattice_vec_a, edge_len, mesh_len, Omega, lattice_type, lattice_constant, Pmaterial ] = ...
%     FAME_Parameter_Generator( grid_nums, material.lattice_type, material.lattice_constant, Pmaterial );
[ Par_mesh, Par_lattice, Par_recip_lattice, Par_material, Par_eig ] = FAME_Parameter_Generator( grid_num, material.lattice_type, material.lattice_constant, Pmaterial);
edge_len      = Par_mesh.edge_len;
mesh_len      = Par_mesh.mesh_len;
Omega         = Par_lattice.Omega;
lattice_type  = Par_lattice.lattice_type;
lattice_constant = Par_lattice.lattice_constant;
Pmaterial = Par_material;

lattice_vec_a = Par_lattice.Omega\Par_lattice.lattice_vec_a;

%% Plot material shape
cla(h.ax_show3D_prim)
cla(h.ax_show3D_unit)

FAME_Main_User_Defined_gui(file_name,Pmaterial.sphere_radius,Pmaterial.cylinder_radius, ...
    lattice_type, lattice_constant, lattice_vec_a, h.ax_show3D_prim, h.ax_show3D_unit);

%% Update Information figure
lattice_constant.a     = num2str(lattice_constant.a);
lattice_constant.b     = num2str(lattice_constant.b);
lattice_constant.c     = num2str(lattice_constant.c);
lattice_constant.alpha = num2str(lattice_constant.alpha);
lattice_constant.beta  = num2str(lattice_constant.beta);
lattice_constant.gamma = num2str(lattice_constant.gamma);
lattice_constant.Permutation = ['[', num2str(lattice_constant.Permutation(1)), ',' , num2str(lattice_constant.Permutation(2)), ',' , num2str(lattice_constant.Permutation(3)),']'];
lattice_vec = {num2str(lattice_vec_a(:,1)'), num2str(lattice_vec_a(:,2)'), num2str(lattice_vec_a(:,3)')};
edge_len    = {num2str(edge_len(1)), num2str(edge_len(2)), num2str(edge_len(3))};
mesh_len    = {num2str(mesh_len(1)), num2str(mesh_len(2)), num2str(mesh_len(3))};
grid_num    = {num2str(grid_num(1)), num2str(grid_num(2)), num2str(grid_num(3))};
sphere_radius   = num2str(Pmaterial.sphere_radius);
cylinder_radius = num2str(Pmaterial.cylinder_radius);

FAME_GUI_Info_Set(  'lattice_vector', lattice_vec, h);
FAME_GUI_Info_Set(     'edge_length', edge_len, h);
FAME_GUI_Info_Set(     'mesh_length', mesh_len, h);
FAME_GUI_Info_Set(    'lattice_type', lattice_type, h);
FAME_GUI_Info_Set('lattice_constant', lattice_constant, h);
FAME_GUI_Info_Set(       'file_name', file_name,  h);
FAME_GUI_Info_Set(    'grid_numbers', grid_num, h);
FAME_GUI_Info_Set(   'sphere_radius', sphere_radius, h);
FAME_GUI_Info_Set( 'cylinder_radius', cylinder_radius, h);
FAME_GUI_Info_Set(    'permittivity', permittivity, h);
FAME_GUI_Info_Set(    'permeability', permeability, h);
FAME_GUI_Info_Set(     'reciprocity', reciprocity, h);
FAME_GUI_Info_Set(       'chirality', chirality, h);
FAME_GUI_Info_Set(   'material_type', material_type, h);
FAME_GUI_Info_Set(        'part_num', part_num, h);

%% Update main FAME figure
set( h.edit_domain.lattice_constant(1), 'string', lattice_constant.a );
set( h.edit_domain.lattice_constant(2), 'string', lattice_constant.b );
set( h.edit_domain.lattice_constant(3), 'string', lattice_constant.c );
set( h.edit_domain.lattice_constant(4), 'string', lattice_constant.alpha );
set( h.edit_domain.lattice_constant(5), 'string', lattice_constant.beta );
set( h.edit_domain.lattice_constant(6), 'string', lattice_constant.gamma );

set( h.edit_domain.grid_num(1), 'string', grid_num{1} );
set( h.edit_domain.grid_num(2), 'string', grid_num{2} );
set( h.edit_domain.grid_num(3), 'string', grid_num{3} );

set( h.edit_domain.shape(1), 'string', sphere_radius );
set( h.edit_domain.shape(2), 'string', cylinder_radius );

set( h.edit_material(1), 'string', permittivity );
set( h.edit_material(2), 'string', permeability );
set( h.edit_material(3), 'string', reciprocity );
set( h.edit_material(4), 'string', chirality );

set( h.edit_lattice_part_num, 'string', part_num );


switch material_type
    case 'isotropic'
        set(h.radio_mode(1), 'value', 1 );
    case 'biisotropic'
        set(h.radio_mode(2), 'value', 1 );
    case 'anisotropic'
        set(h.radio_mode(3), 'value', 1 );
end
