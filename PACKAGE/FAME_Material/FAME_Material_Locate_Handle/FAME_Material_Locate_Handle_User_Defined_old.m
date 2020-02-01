function [ Point_idx ] = FAME_Material_Locate_Handle_User_Defined_old( X,Y,Z,a1,a2,a3,material )

for i = 1:length(material.parameters)
    Point_idx{i} = [];
    % find points inside spheres
    if isfield(material.parameters{i},'sphere_centers') == 1
        for j = 1:size(material.parameters{i}.sphere_centers,1)
            center = ([a1,a2,a3]*material.parameters{i}.sphere_centers(j,:)')';
            a1 = a1';
            a2 = a2';
            a3 = a3';
            center_periodic = [       center;      center+a1;     center+a1+a2;    center+a1+a3; center+a1+a2+a3; ...
                                center+a1-a2;   center+a1-a3;  center+a1-a2-a3; center+a1-a2+a3; center+a1+a2-a3; ...
                                   center-a1;   center-a1+a2;     center-a1+a3; center-a1+a2+a3; ...
                                center-a1-a2;   center-a1-a3;  center-a1-a2-a3; center-a1-a2+a3; center-a1+a2-a3; ...
                                   center+a2;   center+a2+a3;     center+a2-a3;       center-a2; center-a2+a3; center-a2-a3; ...
                                   center+a3;      center-a3] ;
            for k=1:27
                 Point_idx{i} = union(Point_idx{i},FAME_Material_Locate_Sphere( [X,Y,Z], material.parameters{i}.sphere_radius(j), center_periodic(k,:) ));
            end
            a1 = a1';
            a2 = a2';
            a3 = a3';
        end
    end

    % find points inside cylinders
    if isfield(material.parameters{i},'cylinder_bot_centers') == 1
        for j = 1:size(material.parameters{i}.cylinder_bot_centers,1)
            bot_center   = ([a1,a2,a3]*material.parameters{i}.cylinder_bot_centers(j,:)')';
            top_center   = ([a1,a2,a3]*material.parameters{i}.cylinder_top_centers(j,:)')';
            a1 = a1';
            a2 = a2';
            a3 = a3';
            bot_center_periodic = [       bot_center;      bot_center+a1;     bot_center+a1+a2;    bot_center+a1+a3; bot_center+a1+a2+a3; ...
                                    bot_center+a1-a2;   bot_center+a1-a3;  bot_center+a1-a2-a3; bot_center+a1-a2+a3; bot_center+a1+a2-a3; ...
                                       bot_center-a1;   bot_center-a1+a2;     bot_center-a1+a3; bot_center-a1+a2+a3; ...
                                    bot_center-a1-a2;   bot_center-a1-a3;  bot_center-a1-a2-a3; bot_center-a1-a2+a3; bot_center-a1+a2-a3; ...
                                       bot_center+a2;   bot_center+a2+a3;     bot_center+a2-a3;       bot_center-a2; bot_center-a2+a3; bot_center-a2-a3; ...
                                       bot_center+a3;      bot_center-a3] ;
            top_center_periodic = [       top_center;      top_center+a1;     top_center+a1+a2;    top_center+a1+a3; top_center+a1+a2+a3; ...
                                    top_center+a1-a2;   top_center+a1-a3;  top_center+a1-a2-a3; top_center+a1-a2+a3; top_center+a1+a2-a3; ...
                                       top_center-a1;   top_center-a1+a2;     top_center-a1+a3; top_center-a1+a2+a3; ...
                                    top_center-a1-a2;   top_center-a1-a3;  top_center-a1-a2-a3; top_center-a1-a2+a3; top_center-a1+a2-a3; ...
                                       top_center+a2;   top_center+a2+a3;     top_center+a2-a3;       top_center-a2;    top_center-a2+a3; top_center-a2-a3; ...
                                       top_center+a3;      top_center-a3] ;
              for k = 1:27                 
                 Point_idx{i} =  union(Point_idx{i},FAME_Material_Locate_Cylinder( [X,Y,Z], material.parameters{i}.cylinder_radius(j), bot_center_periodic(k,:), top_center_periodic(k,:) ));
              end
              a1 = a1';
            a2 = a2';
            a3 = a3';
        end
    end
end