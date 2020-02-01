function FAME_mfile2txt_future_version( mfile_name, txt_name_para, txt_name_B, txt_name_wave_vec_array )
    % load indicate m-file
    load(mfile_name)
%% Generate text file for B matrix    
    fileID = fopen([txt_name_B,'.txt'],'w');
    fprintf(fileID,'%s \n','# Material inner index (Edge centers for E_1)');
    fprintf(fileID,'%d ',Par.material.B.ele_x_idx{1});
    fprintf(fileID,'\n');
    fprintf(fileID,'%s \n','# Material inner index (Edge centers for E_2)');
    fprintf(fileID,'%d ',Par.material.B.ele_y_idx{1});
    fprintf(fileID,'\n');
    fprintf(fileID,'%s \n','# Material inner index (Edge centers for E_3)');
    fprintf(fileID,'%d ',Par.material.B.ele_z_idx{1});    
    fclose(fileID); 
%% Generate text file for wave vector array
    fileID = fopen([txt_name_wave_vec_array,'.txt'],'w');
    fprintf(fileID,'%s \n','# Total wave vector number');
    fprintf(fileID,'%d \n',Par.recip_lattice.wave_vec_num);
    fprintf(fileID,'%s \n','# wave vector array(a vector of three)');
    fprintf(fileID,'%f ',Par.recip_lattice.wave_vec_array);
    fclose(fileID); 
%% Generate text file for parameters
    fid=fopen(txt_name_para,'wt');
    fprintf(fid,'%s \n','########## Mesh information ##########');
    fprintf(fid,'%s \n','# Grid number');
    fprintf(fid,'%s \n',num2str(Par.mesh.grid_num));
    fprintf(fid,'%s \n','# Edge length');
    fprintf(fid,'%s \n',num2str(Par.mesh.edge_len));
    fprintf(fid,'%s \n','# Mesh length');
    fprintf(fid,'%s \n',num2str(Par.mesh.mesh_len));
    fprintf(fid,'%s \n','########## Lattice information ##########');
    fprintf(fid,'%s \n','# Lattice type');
    fprintf(fid,'%s \n',Par.lattice.lattice_type);
    fprintf(fid,'%s \n','# Lattice constant');
    fprintf(fid,'%s \n',num2str([Par.lattice.lattice_constant.a,...
                                 Par.lattice.lattice_constant.b,...
                                 Par.lattice.lattice_constant.c,...
                                 Par.lattice.lattice_constant.alpha,...
                                 Par.lattice.lattice_constant.beta,...
                                 Par.lattice.lattice_constant.gamma]));
    fprintf(fid,'%s \n','# Lattice vectors(storage in [a_1, a_2, a3])');
    fprintf(fid,'%s \n',num2str(Par.lattice.lattice_vec_a(:)'));
    fprintf(fid,'%s \n','# Rotation matrix(storage in [w_11, w_21, w31, w_12, w_22, w32, w_13, w_23, w33])');
    fprintf(fid,'%s \n',num2str(Par.lattice.Omega(:)'));
    fprintf(fid,'%s \n','########## Reciprocal lattice information ##########');
    fprintf(fid,'%s \n','# Wave vector number');
    fprintf(fid,'%s \n',num2str(Par.recip_lattice.wave_vec_num));
    fprintf(fid,'%s \n','# Brillouin zone Path string');
    fprintf(fid,'%s \n',Par.recip_lattice.path_string);
    fprintf(fid,'%s \n','# Brillouin zone vertex');
    try
        fprintf(fid,'%s ',cell2mat(fieldnames(Par.recip_lattice.vertex)'));
    end
    fprintf(fid,'%s \n','########## Material information ##########');
    fprintf(fid,'%s \n','# Sphere radius');
    if isfield(Par.material,'sphere_radius')
        fprintf(fid,'%s \n',num2str(Par.material.sphere_radius));
    else
        fprintf(fid,'\n');
    end
    fprintf(fid,'%s \n','# Cylinder radius');
    if isfield(Par.material,'cylinder_radius')
        fprintf(fid,'%s \n',num2str(Par.material.cylinder_radius));
    else
        fprintf(fid,'\n');
    end
    fprintf(fid,'%s \n','# Material data name');
    if isfield(Par.material,'data_name')
        fprintf(fid,'%s \n',Par.material.data_name);
    else
        fprintf(fid,'\n');
    end
    fprintf(fid,'%s \n','# Material type');
    fprintf(fid,'%s \n',Par.material.material_type);
    fprintf(fid,'%s \n','# Permittivity(inner material)');
    fprintf(fid,'%s \n',num2str(Par.material.ele_permitt_in));
    fprintf(fid,'%s \n','# Permittivity(outer material)');
    fprintf(fid,'%s \n',num2str(Par.material.ele_permitt_out));
    fprintf(fid,'%s \n','# Permeability(inner material)');
    fprintf(fid,'%s \n',num2str(Par.material.mag_permeab_in));
    fprintf(fid,'%s \n','# Permeability(outer material)');
    fprintf(fid,'%s \n',num2str(Par.material.mag_permeab_out));
    fprintf(fid,'%s \n','# Reciprocity(inner material)');
    if isfield(Par.material,'reciprocity_in')
        fprintf(fid,'%s \n',num2str(Par.material.reciprocity_in));
    else
        fprintf(fid,'\n');
    end
    fprintf(fid,'%s \n','# Reciprocity(outer material)');
    if isfield(Par.material,'reciprocity_out')
        fprintf(fid,'%s \n',num2str(Par.material.reciprocity_out));
    else
        fprintf(fid,'\n');
    end
    fprintf(fid,'%s \n','# Chirality(inner material)');
    if isfield(Par.material,'chirality_in')
        fprintf(fid,'%s \n',num2str(Par.material.chirality_in));
    else
        fprintf(fid,'\n');
    end
    fprintf(fid,'%s \n','# Chirality(outer material)');
    if isfield(Par.material,'chirality_out')
        fprintf(fid,'%s \n',num2str(Par.material.chirality_out));
    else
        fprintf(fid,'\n');
    end
    fprintf(fid,'%s \n','########## Eigensolver information ##########');
    fprintf(fid,'%s \n','# Desired eigenpair number');
    fprintf(fid,'%s \n',num2str(Par.eig.eigen_wanted));
    
    fclose(fid);
end