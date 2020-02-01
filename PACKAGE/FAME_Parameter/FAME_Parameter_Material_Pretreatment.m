function material = FAME_Parameter_Material_Pretreatment(material,pretreatment_idx)
    for i = 1:length(material.parameters)
        if isfield(material.parameters{i},'sphere_radius') == 1
            material.parameters{i}.sphere_centers = material.parameters{i}.sphere_centers(:,pretreatment_idx);
        end
        if isfield(material.parameters{i},'cylinder_radius') == 1
            material.parameters{i}.cylinder_bot_centers = material.parameters{i}.cylinder_bot_centers(:,pretreatment_idx);
            material.parameters{i}.cylinder_top_centers = material.parameters{i}.cylinder_top_centers(:,pretreatment_idx);
        end
    end
end