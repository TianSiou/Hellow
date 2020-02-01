function FAME_Main_User_Defined_gui(data_name,sphere_radius,cylinder_radius, lattice_type, lattice_constant, lattice_vec_a, hax_prim, hax_unit)
%% Select material parameter data
material = [];
eval(['material = ',data_name,'(sphere_radius,cylinder_radius);']);

if isfield(lattice_constant,'Permutation') == 1
    material = FAME_Parameter_Material_Pretreatment(material, lattice_constant.Permutation);
end
%% Plot material shape in primitive cell
FAME_Plot_Material_User_Defined(lattice_vec_a(:,1),lattice_vec_a(:,2),lattice_vec_a(:,3), lattice_type, lattice_constant, material.parameters, hax_prim, 'primitive_cell');
FAME_Plot_Material_User_Defined(lattice_vec_a(:,1),lattice_vec_a(:,2),lattice_vec_a(:,3), lattice_type, lattice_constant, material.parameters, hax_unit, 'unit_cell');
end
