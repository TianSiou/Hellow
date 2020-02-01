function [sphere_radius, cylinder_radius] = FAME_Material_Locate_Parameter_Foolproof(sphere_radius, cylinder_radius, material_num)
    if length(sphere_radius) < material_num
        sphere_radius = sphere_radius*ones(material_num,1);
    end
    if length(cylinder_radius) < material_num
        cylinder_radius = cylinder_radius*ones(material_num,1);
    end
end