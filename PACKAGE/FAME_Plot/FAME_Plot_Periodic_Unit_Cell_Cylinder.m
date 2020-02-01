function [cylinder_bot_centers,cylinder_top_centers,cylinder_radius] = ...
            FAME_Plot_Periodic_Unit_Cell_Cylinder(a1,a2,a3,lattice_type,lattice_constant,cylinder_bot_centers, cylinder_top_centers, cylinder_radius)
    
    switch lattice_type
        case {'simple_cubic','triclinic','primitive_monoclinic'}
            % Not thing need to de in this case.
        case {'face_centered_cubic','face_centered_orthorhombic'}
            unit_cell_edge = sqrt(2)*norm(a1) + 1e-3;

            shift = [ -1  0  0;
                      -1  1  0;
                      -1  0  1;
                       0 -1  0;
                       0 -1  1;
                       0  0 -1;
                       0  0  0;
                       0  0  1;
                       0  1 -1;
                       0  1  0;
                       1  0  0;
                       1 -1  0;
                       1  0 -1 ];

            n_cylinder = length(cylinder_radius);

            periodic_cylinder_bot_centers = [];
            periodic_cylinder_top_centers = [];
            periodic_cylinder_radius      = [];
            for i = 1:size(shift,1)   
                for j = 1:n_cylinder
                    bot_center_temp = cylinder_bot_centers(j,:) + shift(i,:);
                    top_center_temp = cylinder_top_centers(j,:) + shift(i,:);
                    point = 0.5*[a1,a2,a3]*(bot_center_temp + top_center_temp)';
%                     if min( sum( abs( (kron(ones(n_cylinder,1),bot_center_temp) - cylinder_bot_centers)' ) ) ) > 1e-6  &&  ...
%                        min( sum( abs( (kron(ones(n_cylinder,1),top_center_temp) - cylinder_top_centers)' ) ) ) > 1e-6  &&  ...
%                        max([a1,a2,a3]*point') < unit_cell_edge + 1e-3 && min([a1,a2,a3]*point') > -1e-3
%                     if max([a1,a2,a3]*point') < unit_cell_edge + 1e-3 && min([a1,a2,a3]*point') > -1e-3
%                         periodic_cylinder_bot_centers = [periodic_cylinder_bot_centers; bot_center_temp];
%                         periodic_cylinder_top_centers = [periodic_cylinder_top_centers; top_center_temp];
%                         periodic_cylinder_radius      = [periodic_cylinder_radius, cylinder_radius(j)];
%                     end
                    % for A(-a1+a2+a3)(a1-a2+a3)
                    vec_n = cross(-a1+a2+a3,a1-a2+a3);
                    fun_test_A1  = @(x,y,z) (vec_n(1)*x + vec_n(2)*y + vec_n(3)*z)/norm(vec_n);
                    fun_test_A1p = @(x,y,z) fun_test_A1(x-(a1(1)+a2(1)-a3(1)),y-(a1(2)+a2(2)-a3(2)),z-(a1(3)+a2(3)-a3(3)));

                    d(1:2) = [ fun_test_A1(point(1),point(2),point(3)), fun_test_A1p(point(1),point(2),point(3)) ];
                    % for A(-a1+a2+a3)(a1+a2-a3)
                    vec_n = cross(-a1+a2+a3,a1+a2-a3);
                    fun_test_A2  = @(x,y,z) (vec_n(1)*x + vec_n(2)*y + vec_n(3)*z)/norm(vec_n);
                    fun_test_A2p = @(x,y,z) fun_test_A2(x-(a1(1)-a2(1)+a3(1)),y-(a1(2)-a2(2)+a3(2)),z-(a1(3)-a2(3)+a3(3)));

                    d(3:4) = [ fun_test_A2(point(1),point(2),point(3)), fun_test_A2p(point(1),point(2),point(3)) ];
                    % for A(a1-a2+a3)(a1+a2-a3)
                    vec_n = cross(a1-a2+a3,a1+a2-a3);
                    fun_test_A3  = @(x,y,z) (vec_n(1)*x + vec_n(2)*y + vec_n(3)*z)/norm(vec_n);
                    fun_test_A3p = @(x,y,z) fun_test_A3(x-(-a1(1)+a2(1)+a3(1)),y-(-a1(2)+a2(2)+a3(2)),z-(-a1(3)+a2(3)+a3(3)));

                    d(5:6) = [ fun_test_A3(point(1),point(2),point(3)), fun_test_A3p(point(1),point(2),point(3)) ];
%                     if d(1)*d(2) <= 1e-3 && d(3)*d(4) <= 1e-3 && d(5)*d(6) <= 1e-3
%                         periodic_sphere_centers = [periodic_sphere_centers; center_temp];
%                         periodic_sphere_radius  = [periodic_sphere_radius, sphere_radius(j)];
%                     end
                    if d(1)*d(2) <= 1e-3 && d(3)*d(4) <= 1e-3 && d(5)*d(6) <= 1e-3
                        periodic_cylinder_bot_centers = [periodic_cylinder_bot_centers; bot_center_temp];
                        periodic_cylinder_top_centers = [periodic_cylinder_top_centers; top_center_temp];
                        periodic_cylinder_radius  = [periodic_cylinder_radius, cylinder_radius(j)];
                    end
                end
            end

            cylinder_bot_centers = [cylinder_bot_centers;periodic_cylinder_bot_centers];
            cylinder_top_centers = [cylinder_top_centers;periodic_cylinder_top_centers];    
            cylinder_radius      = [cylinder_radius,periodic_cylinder_radius];
        case {'body_centered_cubic','body_centered_tetragonal','body_centered_orthorhombic'}
            unit_cell_edge = (2/sqrt(3))*norm(a1);

            shift = [  0  0  0;
                       1  0  0;
                       0  1  0;
                       0  0  1;
                       1  1  0;
                       1  0  1;
                       0  1  1;
                       1  1  1;];

            n_cylinder     = length(cylinder_radius);
            
            periodic_cylinder_bot_centers = [];
            periodic_cylinder_top_centers = [];
            periodic_cylinder_radius      = [];
            for i = 1:size(shift,1)   
                for j = 1:n_cylinder
                    bot_center_temp = cylinder_bot_centers(j,:) + shift(i,:);
                    top_center_temp = cylinder_top_centers(j,:) + shift(i,:);
                    point = .5*( bot_center_temp + top_center_temp );
                    point = [a1,a2,a3]*point';
                    % for A(a1+a3)(a1+a2)
                    vec_n = cross(a1+a3,a1+a2);
                    fun_test_A1  = @(x,y,z) (vec_n(1)*x + vec_n(2)*y + vec_n(3)*z)/norm(vec_n);
                    fun_test_A1p = @(x,y,z) fun_test_A1(x-a2(1)-a3(1),y-a2(2)-a3(2),z-a2(3)-a3(3));

                    d(1:2) = [ fun_test_A1(point(1),point(2),point(3)), fun_test_A1p(point(1),point(2),point(3)) ];
                    % for A(a2+a3)(a1+a2)
                    vec_n = cross(a2+a3,a1+a2);
                    fun_test_A2  = @(x,y,z) (vec_n(1)*x + vec_n(2)*y + vec_n(3)*z)/norm(vec_n);
                    fun_test_A2p = @(x,y,z) fun_test_A2(x-a1(1)-a3(1),y-a1(2)-a3(2),z-a1(3)-a3(3));

                    d(3:4) = [ fun_test_A2(point(1),point(2),point(3)), fun_test_A2p(point(1),point(2),point(3)) ];
                    % for A(a2+a3)(a1,a3)
                    vec_n = cross(a2+a3,a1+a3);
                    fun_test_A3  = @(x,y,z) (vec_n(1)*x + vec_n(2)*y + vec_n(3)*z)/norm(vec_n);
                    fun_test_A3p = @(x,y,z) fun_test_A3(x-a1(1)-a2(1),y-a1(2)-a2(2),z-a1(3)-a2(3));

                    d(5:6) = [ fun_test_A3(point(1),point(2),point(3)), fun_test_A3p(point(1),point(2),point(3)) ];
                    if d(1)*d(2) <= 1e-3 && d(3)*d(4) <= 1e-3 && d(5)*d(6) <= 1e-3
                        periodic_cylinder_bot_centers = [periodic_cylinder_bot_centers; bot_center_temp];
                        periodic_cylinder_top_centers = [periodic_cylinder_top_centers; top_center_temp];
                        periodic_cylinder_radius      = [periodic_cylinder_radius, cylinder_radius(j)];
                    end
                end
            end
            
            idx_delete     = [];
            for i = 1:n_cylinder
                bot_center_temp = cylinder_bot_centers(i,:);
                top_center_temp = cylinder_top_centers(i,:);
                point = .5*( bot_center_temp + top_center_temp );
                point       = [a1,a2,a3]*point';
%                 if max([a1,a2,a3]*bot_center_temp') > unit_cell_edge + (1e-6) || min([a1,a2,a3]*bot_center_temp') < -1e-6 || ...
%                    max([a1,a2,a3]*top_center_temp') > unit_cell_edge + (1e-6) || min([a1,a2,a3]*top_center_temp') < -1e-6
%                     idx_delete = [idx_delete,i];
%                 end
                d(1:2) = [ fun_test_A1(point(1),point(2),point(3)), fun_test_A1p(point(1),point(2),point(3)) ];
                d(3:4) = [ fun_test_A2(point(1),point(2),point(3)), fun_test_A2p(point(1),point(2),point(3)) ];
                d(5:6) = [ fun_test_A3(point(1),point(2),point(3)), fun_test_A3p(point(1),point(2),point(3)) ];
%                 if d(1)*d(2) >= 1e-3 || d(3)*d(4) >= 1e-3 || d(5)*d(6) >= 1e-3
%                     idx_delete = [idx_delete,i];
%                 end
                if d(1)*d(2) >= 1e-3 || d(3)*d(4) >= 1e-3 || d(5)*d(6) >= 1e-3
                    idx_delete = [idx_delete,i];
                end
%                 if max([a1,a2,a3]*point') > unit_cell_edge + (1e-6) || min([a1,a2,a3]*point') < -1e-6
%                     idx_delete = [idx_delete,i];
%                 end
            end
            cylinder_bot_centers(idx_delete,:) = [];
            cylinder_top_centers(idx_delete,:) = [];
            cylinder_radius(idx_delete)        = [];

            cylinder_bot_centers = [cylinder_bot_centers;periodic_cylinder_bot_centers];
            cylinder_top_centers = [cylinder_top_centers;periodic_cylinder_top_centers];    
            cylinder_radius      = [cylinder_radius,periodic_cylinder_radius];
        case 'hexagonal'
            shift = [  1  0  0;
                       0  1  0;
                       1  1  0];
            shift = shift(:,lattice_constant.Permutation);
            temp_a      = [a1,a2,a3];
            temp_a      = temp_a(:,lattice_constant.Permutation);
            a1_orig     = temp_a(:,1);
            a2_orig     = temp_a(:,2);
            a3_orig     = temp_a(:,3);
            
            n_cylinder = length(cylinder_radius);
            periodic_cylinder_bot_centers       = [];
            periodic_cylinder_top_centers       = [];
            periodic_cylinder_radius            = [];
            
            for i = 1:size(shift,1)   
                for j = 1:n_cylinder
                    bot_center_temp = cylinder_bot_centers(j,:) + shift(i,:);
                    top_center_temp = cylinder_top_centers(j,:) + shift(i,:);
                    point = [a1,a2,a3]*.5*(bot_center_temp+top_center_temp)';

                    % for A13
                    vec_n = cross(a1_orig,a3_orig);
                    fun_test_A13  = @(x,y,z) vec_n(1)*x + vec_n(2)*y + vec_n(3)*z;
                    fun_test_A13p = @(x,y,z) fun_test_A13(x-(a1_orig(1)+2*a2_orig(1)),y-(a1_orig(2)+2*a2_orig(2)),z-(a1_orig(3)+2*a2_orig(3)));
                    d(1:2) = [ fun_test_A13(point(1),point(2),point(3)), fun_test_A13p(point(1),point(2),point(3)) ];
                    % for A23
                    vec_n = cross(a2_orig,a3_orig);
                    fun_test_A23  = @(x,y,z) vec_n(1)*x + vec_n(2)*y + vec_n(3)*z;
                    fun_test_A23p = @(x,y,z) fun_test_A23(x-2*a1_orig(1),y-2*a1_orig(2),z-2*a1_orig(3));
                    d(3:4) = [ fun_test_A23(point(1),point(2),point(3)), fun_test_A23p(point(1),point(2),point(3)) ];
                    % for AA13
                    vec_n = cross(a1_orig + a2_orig,a3_orig);
                    fun_test_AA13  = @(x,y,z) vec_n(1)*(x-a1_orig(1)) + vec_n(2)*(y-a1_orig(2)) + vec_n(3)*(z-a1_orig(3));
                    fun_test_AA13p = @(x,y,z) vec_n(1)*(x-a2_orig(1)) + vec_n(2)*(y-a2_orig(2)) + vec_n(3)*(z-a2_orig(3));
                    d(5:6) = [ fun_test_AA13(point(1),point(2),point(3)), fun_test_AA13p(point(1),point(2),point(3)) ];
                    if d(1)*d(2) < 1e-3 && d(3)*d(4) < 1e-3 && d(5)*d(6) < 1e-3
                        periodic_cylinder_bot_centers = [periodic_cylinder_bot_centers; bot_center_temp];
                        periodic_cylinder_top_centers = [periodic_cylinder_top_centers; top_center_temp];
                        periodic_cylinder_radius      = [periodic_cylinder_radius, cylinder_radius(j)];
                    end
                end
            end
            cylinder_bot_centers = [cylinder_bot_centers;periodic_cylinder_bot_centers];
            cylinder_top_centers = [cylinder_top_centers;periodic_cylinder_top_centers];    
            cylinder_radius      = [cylinder_radius,periodic_cylinder_radius];
        case 'rhombohedral'
            shift = [  0  0 -1;
                       0 -1 -1;
                       0  1 -1;
                       0  1  0;
                       1  1  0;
                       0 -1  0;
                       1 -1  0;
                       1  0  0;
                       1  0 -1];

            n_cylinder = length(cylinder_radius);
            periodic_cylinder_bot_centers       = [];
            periodic_cylinder_top_centers       = [];
            periodic_cylinder_radius            = [];

            for i = 1:size(shift,1)            
                for j = 1:n_cylinder
                    bot_center_temp = cylinder_bot_centers(j,:) + shift(i,:);
                    top_center_temp = cylinder_top_centers(j,:) + shift(i,:);
                    point = [a1,a2,a3]*.5*(bot_center_temp+top_center_temp)';
                    % for A(a1-a2)(a2-a3)
                    vec_n = cross(a1-a2,a2-a3);
                    fun_test_A1  = @(x,y,z) vec_n(1)*x + vec_n(2)*y + vec_n(3)*z;
                    fun_test_A1p = @(x,y,z) fun_test_A1(x-(a1(1)+a2(1)+a3(1)),y-(a1(2)+a2(2)+a3(2)),z-(a1(3)+a2(3)+a3(3)));

                    d(1:2) = [ fun_test_A1(point(1),point(2),point(3)), fun_test_A1p(point(1),point(2),point(3)) ];
                    % for A(a1-a2)(a1+a2+a3)
                    vec_n = cross(a1-a2,a1+a2+a3);
                    fun_test_A2  = @(x,y,z) vec_n(1)*x + vec_n(2)*y + vec_n(3)*z;
                    fun_test_A2p = @(x,y,z) fun_test_A2(x-(a2(1)-a3(1)),y-(a2(2)-a3(2)),z-(a2(3)-a3(3)));

                    d(3:4) = [ fun_test_A2(point(1),point(2),point(3)), fun_test_A2p(point(1),point(2),point(3)) ];
                    % for A(a2-a3)(a1+a2+a3)
                    vec_n = cross(a2-a3,a1+a2+a3);
                    fun_test_A3  = @(x,y,z) vec_n(1)*x + vec_n(2)*y + vec_n(3)*z;
                    fun_test_A3p = @(x,y,z) fun_test_A3(x-(a1(1)-a2(1)),y-(a1(2)-a2(2)),z-(a1(3)-a2(3)));
                    
                    d(5:6) = [ fun_test_A3(point(1),point(2),point(3)), fun_test_A3p(point(1),point(2),point(3)) ];
                    if d(1) >= -1e-3 && d(2) <= 1e-3 && d(3) <= 1e-3 && d(4) >= -1e-3 && d(5) >= -1e-3 && d(6) <= 1e-3
                        periodic_cylinder_bot_centers = [periodic_cylinder_bot_centers; bot_center_temp];
                        periodic_cylinder_top_centers = [periodic_cylinder_top_centers; top_center_temp];
                        periodic_cylinder_radius      = [periodic_cylinder_radius, cylinder_radius(j)];
                    end
                end
            end
            
            idx_delete   = [];
            for i = 1:n_cylinder
                bot_center_temp = cylinder_bot_centers(i,:);
                top_center_temp = cylinder_top_centers(i,:);
                point = [a1,a2,a3]*.5*(bot_center_temp+top_center_temp)';
                d(1:2) = [ fun_test_A1(point(1),point(2),point(3)), fun_test_A1p(point(1),point(2),point(3)) ];
                d(3:4) = [ fun_test_A2(point(1),point(2),point(3)), fun_test_A2p(point(1),point(2),point(3)) ];
                d(5:6) = [ fun_test_A3(point(1),point(2),point(3)), fun_test_A3p(point(1),point(2),point(3)) ];
                if d(1) <= -1e-3 || d(2) >= 1e-3 || d(3) >= 1e-3 || d(4) <= -1e-3 || d(5) <= -1e-3 || d(6) >= 1e-3
                    idx_delete = [idx_delete,i];
                end
            end
            cylinder_bot_centers(idx_delete,:) = [];
            cylinder_top_centers(idx_delete,:) = [];
            cylinder_radius(idx_delete)        = [];
            
            cylinder_bot_centers = [cylinder_bot_centers;periodic_cylinder_bot_centers];
            cylinder_top_centers = [cylinder_top_centers;periodic_cylinder_top_centers];    
            cylinder_radius      = [cylinder_radius,periodic_cylinder_radius];
        case 'a_base_centered_monoclinic'
            shift = [ 0 1 0;
                      0 0 1;
                      0 1 1];
            n_cylinder = length(cylinder_radius);
            periodic_cylinder_bot_centers = [];
            periodic_cylinder_top_centers = [];
            periodic_cylinder_radius      = [];
            for i = 1:size(shift,1)   
                for j = 1:n_cylinder
                    bot_center_temp = cylinder_bot_centers(j,:) + shift(i,:);
                    top_center_temp = cylinder_top_centers(j,:) + shift(i,:);
                    point = [a1,a2,a3]*.5*(bot_center_temp+top_center_temp)';
                    % for A(a3-a2)(a3+a2)
                    vec_n = cross(a3-a2,a3+a2);
                    fun_test_A1  = @(x,y,z) (vec_n(1)*(x - a2(1)) + vec_n(2)*(y - a2(2)) + vec_n(3)*(z - a2(3)))/norm(vec_n);
                    fun_test_A1p = @(x,y,z) (vec_n(1)*(x - a2(1) - a1(1)) + vec_n(2)*(y - a2(2) - a1(2)) + vec_n(3)*(z - a2(3) - a1(3)))/norm(vec_n);

                    d(1:2) = [ fun_test_A1(point(1),point(2),point(3)), fun_test_A1p(point(1),point(2),point(3)) ];
                    % for A(a3-a2)(a1)
                    vec_n = cross(a3-a2,a1);
                    fun_test_A2  = @(x,y,z) (vec_n(1)*(x - a2(1)) + vec_n(2)*(y - a2(2)) + vec_n(3)*(z - a2(3)))/norm(vec_n);
                    fun_test_A2p = @(x,y,z) (vec_n(1)*(x - a3(1) - 2*a2(1)) + vec_n(2)*(y - a3(2) - 2*a2(2)) + vec_n(3)*(z - a3(3) - 2*a2(3)))/norm(vec_n);

                    d(3:4) = [ fun_test_A2(point(1),point(2),point(3)), fun_test_A2p(point(1),point(2),point(3)) ];
                    % for A(a3+a2)(a1)
                    vec_n = cross(a3+a2,a1);
                    fun_test_A3  = @(x,y,z) (vec_n(1)*(x - a2(1)) + vec_n(2)*(y - a2(2)) + vec_n(3)*(z - a2(3)))/norm(vec_n);
                    fun_test_A3p = @(x,y,z) (vec_n(1)*(x - a3(1)) + vec_n(2)*(y - a3(2)) + vec_n(3)*(z - a3(3)))/norm(vec_n);

                    d(5:6) = [ fun_test_A3(point(1),point(2),point(3)), fun_test_A3p(point(1),point(2),point(3)) ];
%                     if d(1)*d(2) <= 1e-3 && d(3)*d(4) <= 1e-3 && d(5)*d(6) <= 1e-3
%                         periodic_cylinder_bot_centers = [periodic_cylinder_bot_centers; bot_center_temp];
%                         periodic_cylinder_top_centers = [periodic_cylinder_top_centers; top_center_temp];
%                         periodic_cylinder_radius      = [periodic_cylinder_radius, cylinder_radius(j)];
%                     end
                    if d(1)*d(2) <= .5 && d(3)*d(4) <= .5 && d(5)*d(6) <= .5
                        periodic_cylinder_bot_centers = [periodic_cylinder_bot_centers; bot_center_temp];
                        periodic_cylinder_top_centers = [periodic_cylinder_top_centers; top_center_temp];
                        periodic_cylinder_radius      = [periodic_cylinder_radius, cylinder_radius(j)];
                    end
                end
            end
            
            idx_delete   = [];
            for i = 1:n_cylinder
                bot_center_temp = cylinder_bot_centers(i,:);
                top_center_temp = cylinder_top_centers(i,:);
                point = [a1,a2,a3]*.5*(bot_center_temp+top_center_temp)';
                d(1:2) = [ fun_test_A1(point(1),point(2),point(3)), fun_test_A1p(point(1),point(2),point(3)) ];
                d(3:4) = [ fun_test_A2(point(1),point(2),point(3)), fun_test_A2p(point(1),point(2),point(3)) ];
                d(5:6) = [ fun_test_A3(point(1),point(2),point(3)), fun_test_A3p(point(1),point(2),point(3)) ];
%                 if d(1)*d(2) >= 1e-3 || d(3)*d(4) >= 1e-3 || d(5)*d(6) >= 1e-3
%                     idx_delete = [idx_delete,i];
%                 end
                if d(1)*d(2) >= 1 || d(3)*d(4) >= 1 || d(5)*d(6) >= 1
                    idx_delete = [idx_delete,i];
                end
            end
            cylinder_bot_centers(idx_delete,:) = [];
            cylinder_top_centers(idx_delete,:) = [];
            cylinder_radius(idx_delete)        = [];
            
            cylinder_bot_centers = [cylinder_bot_centers;periodic_cylinder_bot_centers];
            cylinder_top_centers = [cylinder_top_centers;periodic_cylinder_top_centers];    
            cylinder_radius      = [cylinder_radius,periodic_cylinder_radius];
        case 'a_base_centered_orthorhombic'
            shift = [ 0 -1 0;
                      0  0 1;
                      0 -1 1];
            n_cylinder = length(cylinder_radius);
            periodic_cylinder_bot_centers = [];
            periodic_cylinder_top_centers = [];
            periodic_cylinder_radius      = [];
            for i = 1:size(shift,1)   
                for j = 1:n_cylinder
                    bot_center_temp = cylinder_bot_centers(j,:) + shift(i,:);
                    top_center_temp = cylinder_top_centers(j,:) + shift(i,:);
                    point = [a1,a2,a3]*.5*(bot_center_temp+top_center_temp)';
                    % for A(a1)(a3+a2)
                    vec_n = cross(a1,a3+a2);
                    fun_test_A1  = @(x,y,z) (vec_n(1)*x + vec_n(2)*y + vec_n(3)*z)/norm(vec_n);
                    fun_test_A1p = @(x,y,z) (vec_n(1)*(x - a3(1) + a2(1)) + vec_n(2)*(y - a3(2) + a2(2)) + vec_n(3)*(z - a3(3) + a2(3)))/norm(vec_n);

                    d(1:2) = [ fun_test_A1(point(1),point(2),point(3)), fun_test_A1p(point(1),point(2),point(3)) ];
                    % for A(a3-a2)(a1)
                    vec_n = cross(a3-a2,a1);
                    fun_test_A2  = @(x,y,z) (vec_n(1)*x + vec_n(2)*y + vec_n(3)*z)/norm(vec_n);
                    fun_test_A2p = @(x,y,z) (vec_n(1)*(x - a3(1) - a2(1)) + vec_n(2)*(y - a3(2) - a2(2)) + vec_n(3)*(z - a3(3) - a2(3)))/norm(vec_n);

                    d(3:4) = [ fun_test_A2(point(1),point(2),point(3)), fun_test_A2p(point(1),point(2),point(3)) ];
                    % for A(a3+a2)(a3-a2)
                    vec_n = cross(a3+a2,a3-a2);
                    fun_test_A3  = @(x,y,z) (vec_n(1)*x + vec_n(2)*y + vec_n(3)*z)/norm(vec_n);
                    fun_test_A3p = @(x,y,z) (vec_n(1)*(x - a1(1)) + vec_n(2)*(y - a1(2)) + vec_n(3)*(z - a1(3)))/norm(vec_n);

                    d(5:6) = [ fun_test_A3(point(1),point(2),point(3)), fun_test_A3p(point(1),point(2),point(3)) ];
                    if d(1)*d(2) <= .5 && d(3)*d(4) <= .5 && d(5)*d(6) <= .5
                        periodic_cylinder_bot_centers = [periodic_cylinder_bot_centers; bot_center_temp];
                        periodic_cylinder_top_centers = [periodic_cylinder_top_centers; top_center_temp];
                        periodic_cylinder_radius      = [periodic_cylinder_radius, cylinder_radius(j)];
                    end
                end
            end
            
            idx_delete   = [];
            for i = 1:n_cylinder
                bot_center_temp = cylinder_bot_centers(i,:);
                top_center_temp = cylinder_top_centers(i,:);
                point = [a1,a2,a3]*.5*(bot_center_temp+top_center_temp)';
                d(1:2) = [ fun_test_A1(point(1),point(2),point(3)), fun_test_A1p(point(1),point(2),point(3)) ];
                d(3:4) = [ fun_test_A2(point(1),point(2),point(3)), fun_test_A2p(point(1),point(2),point(3)) ];
                d(5:6) = [ fun_test_A3(point(1),point(2),point(3)), fun_test_A3p(point(1),point(2),point(3)) ];
%                 if d(1)*d(2) >= 1e-3 || d(3)*d(4) >= 1e-3 || d(5)*d(6) >= 1e-3
%                     idx_delete = [idx_delete,i];
%                 end
                if d(1)*d(2) >= 1 || d(3)*d(4) >= 1 || d(5)*d(6) >= 1
                    idx_delete = [idx_delete,i];
                end
            end
            cylinder_bot_centers(idx_delete,:) = [];
            cylinder_top_centers(idx_delete,:) = [];
            cylinder_radius(idx_delete)        = [];
            
            cylinder_bot_centers = [cylinder_bot_centers;periodic_cylinder_bot_centers];
            cylinder_top_centers = [cylinder_top_centers;periodic_cylinder_top_centers];    
            cylinder_radius      = [cylinder_radius,periodic_cylinder_radius];
        case 'c_base_centered_orthorhombic'
            shift = [  0  1 0;
                      -1  0 0;
                      -1  1 0];
            n_cylinder = length(cylinder_radius);
            periodic_cylinder_bot_centers = [];
            periodic_cylinder_top_centers = [];
            periodic_cylinder_radius      = [];
            for i = 1:size(shift,1)   
                for j = 1:n_cylinder
                    bot_center_temp = cylinder_bot_centers(j,:) + shift(i,:);
                    top_center_temp = cylinder_top_centers(j,:) + shift(i,:);
                    point = [a1,a2,a3]*.5*(bot_center_temp+top_center_temp)';
                    % for A(a3)(a1+a2)
                    vec_n = cross(a3,a1+a2);
                    fun_test_A1  = @(x,y,z) (vec_n(1)*x + vec_n(2)*y + vec_n(3)*z)/norm(vec_n);
                    fun_test_A1p = @(x,y,z) (vec_n(1)*(x - a2(1) + a1(1)) + vec_n(2)*(y - a2(2) + a1(2)) + vec_n(3)*(z - a2(3) + a1(3)))/norm(vec_n);

                    d(1:2) = [ fun_test_A1(point(1),point(2),point(3)), fun_test_A1p(point(1),point(2),point(3)) ];
                    % for A(a2-a1)(a3)
                    vec_n = cross(a2-a1,a3);
                    fun_test_A2  = @(x,y,z) (vec_n(1)*x + vec_n(2)*y + vec_n(3)*z)/norm(vec_n);
                    fun_test_A2p = @(x,y,z) (vec_n(1)*(x - a2(1) - a1(1)) + vec_n(2)*(y - a2(2) - a1(2)) + vec_n(3)*(z - a2(3) - a1(3)))/norm(vec_n);

                    d(3:4) = [ fun_test_A2(point(1),point(2),point(3)), fun_test_A2p(point(1),point(2),point(3)) ];
                    % for A(a1+a2)(a2-a1)
                    vec_n = cross(a2+a1,a2-a1);
                    fun_test_A3  = @(x,y,z) (vec_n(1)*x + vec_n(2)*y + vec_n(3)*z)/norm(vec_n);
                    fun_test_A3p = @(x,y,z) (vec_n(1)*(x - a3(1)) + vec_n(2)*(y - a3(2)) + vec_n(3)*(z - a3(3)))/norm(vec_n);
                    
                    d(5:6) = [ fun_test_A3(point(1),point(2),point(3)), fun_test_A3p(point(1),point(2),point(3)) ];
                    if d(1)*d(2) <= 1e-3 && d(3)*d(4) <= 1e-3 && d(5)*d(6) <= 1e-3
                        periodic_cylinder_bot_centers = [periodic_cylinder_bot_centers; bot_center_temp];
                        periodic_cylinder_top_centers = [periodic_cylinder_top_centers; top_center_temp];
                        periodic_cylinder_radius      = [periodic_cylinder_radius, cylinder_radius(j)];
                    end
                end
            end
            
            idx_delete   = [];
            for i = 1:n_cylinder
                bot_center_temp = cylinder_bot_centers(i,:);
                top_center_temp = cylinder_top_centers(i,:);
                point = [a1,a2,a3]*.5*(bot_center_temp+top_center_temp)';
                d(1:2) = [ fun_test_A1(point(1),point(2),point(3)), fun_test_A1p(point(1),point(2),point(3)) ];
                d(3:4) = [ fun_test_A2(point(1),point(2),point(3)), fun_test_A2p(point(1),point(2),point(3)) ];
                d(5:6) = [ fun_test_A3(point(1),point(2),point(3)), fun_test_A3p(point(1),point(2),point(3)) ];
%                 if d(1)*d(2) >= 1e-3 || d(3)*d(4) >= 1e-3 || d(5)*d(6) >= 1e-3
%                     idx_delete = [idx_delete,i];
%                 end
                if d(1)*d(2) >= -1e-3 || d(3)*d(4) >= -1e-3 || d(5)*d(6) >= -1e-3
                    idx_delete = [idx_delete,i];
                end
            end
            cylinder_bot_centers(idx_delete,:) = [];
            cylinder_top_centers(idx_delete,:) = [];
            cylinder_radius(idx_delete)        = [];
            
            cylinder_bot_centers = [cylinder_bot_centers;periodic_cylinder_bot_centers];
            cylinder_top_centers = [cylinder_top_centers;periodic_cylinder_top_centers];    
            cylinder_radius      = [cylinder_radius,periodic_cylinder_radius];
    end
end