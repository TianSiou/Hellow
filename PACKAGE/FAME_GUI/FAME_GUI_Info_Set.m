function FAME_GUI_Info_Set(parameter_name, set_value, h)
    switch parameter_name
        case 'file_name'
            set(h.info_text_material_file_value,'string',set_value);
        case 'material_type'
            set(h.info_text_material_type_value,'string',set_value); 
        case 'grid_numbers'
            set(h.info_text_grid_num_value(1),'string',set_value{1});
            set(h.info_text_grid_num_value(2),'string',set_value{2});
            set(h.info_text_grid_num_value(3),'string',set_value{3});
        case 'mesh_length'
            set(h.info_text_mesh_len_value(1),'string',set_value{1});
            set(h.info_text_mesh_len_value(2),'string',set_value{2});
            set(h.info_text_mesh_len_value(3),'string',set_value{3});
        case 'edge_length'
            set(h.info_text_edge_len_value(1),'string',set_value{1});
            set(h.info_text_edge_len_value(2),'string',set_value{2});
            set(h.info_text_edge_len_value(3),'string',set_value{3});
        case 'lattice_type'
            set(h.info_text_lattice_type_value,'string',set_value);    
        case 'lattice_vector'
            set(h.info_text_lattice_vector_value(1),'string',set_value{1});
            set(h.info_text_lattice_vector_value(2),'string',set_value{2});
            set(h.info_text_lattice_vector_value(3),'string',set_value{3});    
        case 'lattice_constant'
            set(h.info_text_lattice_constant_value(1),'string',set_value.a);
            set(h.info_text_lattice_constant_value(2),'string',set_value.b);
            set(h.info_text_lattice_constant_value(3),'string',set_value.c);
            set(h.info_text_lattice_constant_value(4),'string',set_value.alpha);
            set(h.info_text_lattice_constant_value(5),'string',set_value.beta);
            set(h.info_text_lattice_constant_value(6),'string',set_value.gamma);
            set(h.info_text_lattice_constant_value(7),'string',set_value.Permutation);
        case 'sphere_radius'
            set(h.info_text_sphere_radius_value,'string',set_value);
        case 'cylinder_radius'
            set(h.info_text_cylinder_radius_value,'string',set_value);    
        case 'permittivity'
            set(h.info_text_material_para_value(1),'string',set_value);    
        case 'permeability'
            set(h.info_text_material_para_value(2),'string',set_value);        
        case 'reciprocity'
            set(h.info_text_material_para_value(3),'string',set_value);        
        case 'chirality'
            set(h.info_text_material_para_value(4),'string',set_value);            
        case 'part_num'    
            set(h.info_text_partition_num_val,'string',set_value);            
    end
end