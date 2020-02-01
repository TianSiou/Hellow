function [Popt] = FAME_GUI_Info_Get(h)
%% Lattice info
    Popt.lattice.lattice_type       = get(h.info_text_lattice_type_value,'string');
    
    Popt.lattice.lattice_constant.a           = str2num( get(h.info_text_lattice_constant_value(1),'string') );
    Popt.lattice.lattice_constant.b           = str2num( get(h.info_text_lattice_constant_value(2),'string') );
    Popt.lattice.lattice_constant.c           = str2num( get(h.info_text_lattice_constant_value(3),'string') );
    Popt.lattice.lattice_constant.alpha       = str2num( get(h.info_text_lattice_constant_value(4),'string') );
    Popt.lattice.lattice_constant.beta        = str2num( get(h.info_text_lattice_constant_value(5),'string') );
    Popt.lattice.lattice_constant.gamma       = str2num( get(h.info_text_lattice_constant_value(6),'string') );
    Popt.lattice.lattice_constant.Permutation = str2num( get(h.info_text_lattice_constant_value(7),'string') );
    
    a1 = reshape( str2num( get(h.info_text_lattice_vector_value(1),'string') ), 3,1);
    a2 = reshape( str2num( get(h.info_text_lattice_vector_value(2),'string') ), 3,1);
    a3 = reshape( str2num( get(h.info_text_lattice_vector_value(3),'string') ), 3,1);
    Popt.lattice.lattice_vector_a = [a1,a2,a3];
    
    Popt.recip_lattice.path_string  = get( h.info_text_BZ_path_value, 'string');
    Popt.recip_lattice.part_num     = str2num(get( h.info_text_partition_num_val, 'string'));
%% Mesh info    
    Popt.mesh.grid_num(1)  = str2num(get( h.info_text_grid_num_value(1),'string'));
    Popt.mesh.grid_num(2)  = str2num(get( h.info_text_grid_num_value(2),'string'));
    Popt.mesh.grid_num(3)  = str2num(get( h.info_text_grid_num_value(3),'string'));
    
    Popt.mesh.edge_len(1)  = str2num(get( h.info_text_edge_len_value(1),'string'));
    Popt.mesh.edge_len(2)  = str2num(get( h.info_text_edge_len_value(2),'string'));
    Popt.mesh.edge_len(3)  = str2num(get( h.info_text_edge_len_value(3),'string'));

    Popt.mesh.mesh_len(1)  = str2num(get( h.info_text_mesh_len_value(1),'string'));
    Popt.mesh.mesh_len(2)  = str2num(get( h.info_text_mesh_len_value(2),'string'));
    Popt.mesh.mesh_len(3)  = str2num(get( h.info_text_mesh_len_value(3),'string'));
%% Material info
    Popt.material.material_type = get( h.info_text_material_type_value,'string');
    
    Popt.material.ele_permitt_in = str2num(get( h.info_text_material_para_value(1),'string'));
    Popt.material.mag_permeab_in = str2num(get( h.info_text_material_para_value(2),'string'));
    Popt.material.reciprocity_in = str2num(get( h.info_text_material_para_value(3),'string'));
    Popt.material.chirality_in   = str2num(get( h.info_text_material_para_value(4),'string'));
    Popt.material.ele_permitt_out = 1;
    Popt.material.mag_permeab_out = 1;
    Popt.material.reciprocity_out = 0;
    Popt.material.chirality_out   = 0;

    Popt.material.data_name   = get( h.info_text_material_file_value,'string');
    Popt.material.sphere_radius   = str2num(get( h.info_text_sphere_radius_value,'string'));
    Popt.material.cylinder_radius = str2num(get( h.info_text_cylinder_radius_value,'string'));
    
    Popt.recip_lattice.part_num = str2num(get( h.info_text_partition_num_val,'string'));
end