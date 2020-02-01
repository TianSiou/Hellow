function [ lattice_type, lattice_constant ] = Data_File_Check(material)

    if isfield(material,'lattice_type') == 0
        error('There are no lattice_type in your data file!');
    else
%         lattice_type = material.lattice_type;
        lattice_type = 'user_defined';
    end

    if isfield(material,'lattice_constant') == 0
        error('There are no lattice_constant in your data file!');
    else
        if isfield(material.lattice_constant,'a')
            lattice_constant.a = material.lattice_constant.a;
        end
        if isfield(material.lattice_constant,'b')
            lattice_constant.b = material.lattice_constant.b;
        end
        if isfield(material.lattice_constant,'c')
            lattice_constant.c = material.lattice_constant.c;
        end
        if isfield(material.lattice_constant,'alpha')
            lattice_constant.alpha = material.lattice_constant.alpha;
        end
        if isfield(material.lattice_constant,'beta')
            lattice_constant.beta = material.lattice_constant.beta;
        end
        if isfield(material.lattice_constant,'gamma')
            lattice_constant.gamma = material.lattice_constant.gamma;
        end
    end
    
    for i = 1:length(material.parameters)
        if isfield(material.parameters{i},'sphere_centers') == 1
            if sum(material.parameters{i}.sphere_centers < 0) > 0
                error(['The sphere centers must be set in [0,1]x[0,1]x[0,1]! There are some error in material.parameters{','num2str(i)','}']);
            end
        end
        if isfield(material.parameters{i},'cylinder_bot_centers') == 1
%             if sum(material.parameters{i}.cylinder_bot_centers < 0) > 0
%                 error(['The cylinder bot centers must be set in [0,1]x[0,1]x[0,1]! There are some error in material.parameters{','num2str(i)','}']);
%             end
%             if sum(material.parameters{i}.cylinder_top_centers < 0) > 0
%                 error(['The cylinder top centers must be set in [0,1]x[0,1]x[0,1]! There are some error in material.parameters{','num2str(i)','}']);
%             end
        end
        
    end

    
end