function [ Point_idx ] = FAME_Material_Locate_Handle_User_Defined( varargin )
    
    X  = varargin{1}; Y  = varargin{2}; Z  = varargin{3};
    a1 = varargin{4}; a2 = varargin{5}; a3 = varargin{6};
    switch nargin
        case 8
            Isofunction     = varargin{7};
            Isovalue        = varargin{8};
            Point_idx       = Locate_by_Isofunction(X,Y,Z,a1,a2,a3,Isofunction,Isovalue);
        case 9
            data_name       = varargin{7};
            sphere_radius   = varargin{8};
            cylinder_radius = varargin{9};
            Point_idx       = Locate_by_SphereCylinder(X,Y,Z,a1,a2,a3,data_name,sphere_radius,cylinder_radius);
        otherwise 
            error('Too many or too few input parameters.');
    end
end
function Point_idx = Locate_by_SphereCylinder(X,Y,Z,a1,a2,a3,data_name,sphere_radius,cylinder_radius)
    Pmaterial = feval(data_name,sphere_radius,cylinder_radius);
    Point_set = [X,Y,Z];
    for i = 1:length(Pmaterial.parameters)
        Point_idx{i} = [];
        % find points inside spheres
        if isfield(Pmaterial.parameters{i},'sphere_centers') == 1
            for j = 1:size(Pmaterial.parameters{i}.sphere_centers,1)
                center = Pmaterial.parameters{i}.sphere_centers(j,:)*[a1,a2,a3]';
                Point_idx{i} = union(Point_idx{i},FAME_Material_Locate_Sphere( Point_set, Pmaterial.parameters{i}.sphere_radius(j), center ));
            end
        end

        % find points inside cylinders
        if isfield(Pmaterial.parameters{i},'cylinder_bot_centers') == 1
            for j = 1:size(Pmaterial.parameters{i}.cylinder_bot_centers,1)
              bot_center = Pmaterial.parameters{i}.cylinder_bot_centers(j,:)*[a1,a2,a3]';
              top_center = Pmaterial.parameters{i}.cylinder_top_centers(j,:)*[a1,a2,a3]';
              Point_idx{i} =  union(Point_idx{i},FAME_Material_Locate_Cylinder( Point_set, Pmaterial.parameters{i}.cylinder_radius(j), bot_center, top_center ));
            end
        end
    end
end
% function Point_idx = Locate_by_SphereCylinder(X,Y,Z,a1,a2,a3,data_name,sphere_radius,cylinder_radius)
%     Pmaterial = feval(data_name,sphere_radius,cylinder_radius);
% 
%     shift = [ -1, -1, -1; 
%                0, -1, -1;
%                1, -1, -1;
%               -1,  0, -1; 
%                0,  0, -1;
%                1,  0, -1;
%               -1,  1, -1; 
%                0,  1, -1;
%                1,  1, -1; 
%               -1, -1,  0; 
%                0, -1,  0;
%                1, -1,  0;
%               -1,  0,  0; 
%                0,  0,  0;
%                1,  0,  0;
%               -1,  1,  0; 
%                0,  1,  0;
%                1,  1,  0;
%               -1, -1,  1; 
%                0, -1,  1;
%                1, -1,  1;
%               -1,  0,  1; 
%                0,  0,  1;
%                1,  0,  1;
%               -1,  1,  1; 
%                0,  1,  1;
%                1,  1,  1];
% 
%     for i = 1:length(Pmaterial.parameters)
%         Point_idx{i} = [];
%         % find points inside spheres
%         if isfield(Pmaterial.parameters{i},'sphere_centers') == 1
%             for j = 1:size(Pmaterial.parameters{i}.sphere_centers,1)
%                 for k=1:size(shift,1)
%                     center_shift = ( Pmaterial.parameters{i}.sphere_centers(j,:) + shift(k,:) )*[a1,a2,a3]';
%                     Point_idx{i} = union(Point_idx{i},FAME_Material_Locate_Sphere( [X,Y,Z], Pmaterial.parameters{i}.sphere_radius(j), center_shift ));
%                 end
%             end
%         end
% 
%         % find points inside cylinders
%         if isfield(Pmaterial.parameters{i},'cylinder_bot_centers') == 1
%             for j = 1:size(Pmaterial.parameters{i}.cylinder_bot_centers,1)
%                 for k = 1:size(shift,1)
%                   bot_center_shift = ( Pmaterial.parameters{i}.cylinder_bot_centers(j,:) + shift(k,:) )*[a1,a2,a3]';
%                   top_center_shift = ( Pmaterial.parameters{i}.cylinder_top_centers(j,:) + shift(k,:) )*[a1,a2,a3]';
%                   Point_idx{i} =  union(Point_idx{i},FAME_Material_Locate_Cylinder( [X,Y,Z], Pmaterial.parameters{i}.cylinder_radius(j), bot_center_shift, top_center_shift ));
%                 end
%             end
%         end
%     end
% end

function Point_idx = Locate_by_Isofunction(X,Y,Z,a1,a2,a3,Isofunction,Isovalue)
    Value = Isofunction(X,Y,Z,a1,a2,a3);
    for i = 1:length(Value)
        Point_idx{i} = find( Value{i} >= Isovalue(i) );
    end
end