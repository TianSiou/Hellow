clear all;
clc;

sphere_radius = 0.1;
cylinder_radius = 0.1;
[ pwd,'\No*.m']
para_file = dir([ pwd,'\No*.m']);
for ii = 1:length(para_file) % ï¿½]ï¿½wï¿½ï¿½ï¿?shape ï¿½Îªï¿½ï¿½Uï¿½Ô¦ï¿½ï¿½ï¿½ï¿?
    dataname = para_file(ii).name(1:end-2)
% dataname  = 'No227_FCC_Cu2Mg';
    eval( ['material  =  ',dataname,'(sphere_radius, cylinder_radius);'] );
    % material  =  No139_BCT_Si(sphere_radius, cylinder_radius);
    fid = fopen([dataname '.txt'], 'w');

    %% Data_name
    fprintf(fid,'#data_name\r\n');
    fprintf(fid, '%s\r\n', dataname);

    %% Material_num
    fprintf(fid,'#material_num\r\n');
    fprintf(fid,'%d\r\n', material.material_num);

    %% Lattice_type
    fprintf(fid,'#lattice_type\r\n');
    fprintf(fid,'%s\r\n', material.lattice_type);

    %% Lattice_constant
    fprintf(fid,'#lattice_constant\r\n');
    if isfield(material.lattice_constant,'a')
        fprintf(fid,'%f ', material.lattice_constant.a);
    else
        fprintf(fid,'0.0 ');
    end
    if isfield(material.lattice_constant,'b')
        fprintf(fid,'%f ', material.lattice_constant.b);
    else
        fprintf(fid,'0.0 ');
    end
    if isfield(material.lattice_constant,'c')
        fprintf(fid,'%f ', material.lattice_constant.c);
    else
        fprintf(fid,'0.0 ');
    end
    if isfield(material.lattice_constant,'alpha')
        fprintf(fid,'%f ', material.lattice_constant.alpha);
    else
        fprintf(fid,'0.0 ');
    end
    if isfield(material.lattice_constant,'beta')
        fprintf(fid,'%f ', material.lattice_constant.beta);
    else
        fprintf(fid,'0.0 ');
    end
    if isfield(material.lattice_constant,'gamma')
        fprintf(fid,'%f ', material.lattice_constant.gamma);
    else
        fprintf(fid,'0.0 ');
    end
    fprintf(fid,'\r\n');

    %% Sphere_num
    fprintf(fid,'#sphere_num\r\n');
    for i = 1 : material.material_num
        if isfield(material.parameters{i},'sphere_radius')
            fprintf(fid, '%d\r\n', size(material.parameters{i}.sphere_radius, 2));
        end
    end

    %% Sphere_centers
    fprintf(fid,'#sphere_centers\r\n');
    for i = 1 : material.material_num
        if isfield(material.parameters{i},'sphere_radius')
            for j = 1 : size(material.parameters{i}.sphere_radius, 2)
                fprintf(fid, '%f %f %f ', material.parameters{i}.sphere_centers(j, :));
            end
            fprintf(fid,'\r\n');
        end
    end

    %% Sphere_radius
    fprintf(fid,'#sphere_radius\r\n');
    for i = 1 : material.material_num
        if isfield(material.parameters{i},'sphere_radius')
            for j = 1 : size(material.parameters{i}.sphere_radius, 2)
                fprintf(fid, '%f ', material.parameters{i}.sphere_radius(j));
            end
            fprintf(fid,'\r\n');
        end
    end


    %% Cylinder_num
    fprintf(fid,'#cylinder_num\r\n');
    for i = 1 : material.material_num
        if isfield(material.parameters{i},'cylinder_radius')
            fprintf(fid, '%d\r\n', size(material.parameters{i}.cylinder_radius, 2));
        else
            fprintf(fid, '0\r\n');
        end
    end

    %% Cylinder_top_centers
    fprintf(fid,'#cylinder_top_centers\r\n');
    for i = 1 : material.material_num
        if isfield(material.parameters{i},'cylinder_radius')
            for j = 1 : size(material.parameters{i}.cylinder_radius, 2)
                fprintf(fid, '%f %f %f ', material.parameters{i}.cylinder_top_centers(j, :));
            end
            fprintf(fid,'\r\n');
        end
    end

    %% Cylinder_bot_centers
    fprintf(fid,'#cylinder_bot_centers\r\n');
    for i = 1 : material.material_num
        if isfield(material.parameters{i},'cylinder_radius')
            for j = 1 : size(material.parameters{i}.cylinder_radius, 2)
                fprintf(fid, '%f %f %f ', material.parameters{i}.cylinder_bot_centers(j, :));
            end
            fprintf(fid,'\r\n');
        end
    end


    %% Cylinder_radius
    fprintf(fid,'#cylinder_radius\r\n');
    for i = 1 : material.material_num
        if isfield(material.parameters{i},'cylinder_radius')
            for j = 1 : size(material.parameters{i}.cylinder_radius, 2)
                fprintf(fid, '%f ', material.parameters{i}.cylinder_radius(j));
            end
            fprintf(fid,'\r\n');
        end
    end

    fclose(fid);
end