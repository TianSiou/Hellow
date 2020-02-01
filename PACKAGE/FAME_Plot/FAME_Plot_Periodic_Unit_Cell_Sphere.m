function [sphere_centers,sphere_radius] = ...
            FAME_Plot_Periodic_Unit_Cell_Sphere(a1,a2,a3,lattice_type,lattice_constant,sphere_centers,sphere_radius)
    switch lattice_type
        case {'simple_cubic','triclinic','primitive_monoclinic','primitive_tetragonal','primitive_othorhombic'}
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
            n_sphere   = length(sphere_radius);

            periodic_sphere_centers       = [];
            periodic_sphere_radius        = [];
            for i = 1:size(shift,1)   
                for j = 1:n_sphere
%                     center_temp = sphere_centers(j,:) + shift(i,:);
%                     if min( sum( abs( (kron(ones(n_sphere,1),center_temp) - sphere_centers)' ) ) ) > 1e-6 && ...
%                        max([a1,a2,a3]*center_temp') < unit_cell_edge && min([a1,a2,a3]*center_temp') > -1e-3;
% 
%                         periodic_sphere_centers = [periodic_sphere_centers; center_temp];
%                         periodic_sphere_radius  = [periodic_sphere_radius, sphere_radius(j)];
%                     end
                    center_temp = sphere_centers(j,:) + shift(i,:);
                    point       = [a1,a2,a3]*center_temp';
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
                        periodic_sphere_centers = [periodic_sphere_centers; center_temp];
                        periodic_sphere_radius  = [periodic_sphere_radius, sphere_radius(j)];
                    end
                end
            end
            sphere_centers       = [sphere_centers;periodic_sphere_centers];
            sphere_radius        = [sphere_radius,periodic_sphere_radius];
        case {'body_centered_cubic','body_centered_tetragonal','body_centered_orthorhombic'}
            unit_cell_edge = (2/sqrt(3))*norm(a1);
            
            shift = [  0  0  0;
                       1  0  0;
                       0  1  0;
                       0  0  1;
                       1  1  0;
                       1  0  1;
                       0  1  1;
                       1  1  1];
                   
            n_sphere = length(sphere_radius);
            periodic_sphere_centers       = [];
            periodic_sphere_radius        = [];
            for i = 1:size(shift,1)   
                for j = 1:n_sphere
                    center_temp = sphere_centers(j,:) + shift(i,:);
                    point       = [a1,a2,a3]*center_temp';
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
%                     if d(1)*d(2) <= 1e-3 && d(3)*d(4) <= 1e-3 && d(5)*d(6) <= 1e-3
%                         periodic_sphere_centers = [periodic_sphere_centers; center_temp];
%                         periodic_sphere_radius  = [periodic_sphere_radius, sphere_radius(j)];
%                     end
                    if d(1)*d(2) <= 1e-3 && d(3)*d(4) <= 1e-3 && d(5)*d(6) <= 1e-3
                        periodic_sphere_centers = [periodic_sphere_centers; center_temp];
                        periodic_sphere_radius  = [periodic_sphere_radius, sphere_radius(j)];
                    end
                end
            end
            
            idx_delete   = [];
            for i = 1:n_sphere
                center_temp = sphere_centers(i,:);
                point       = [a1,a2,a3]*center_temp';
                d(1:2) = [ fun_test_A1(point(1),point(2),point(3)), fun_test_A1p(point(1),point(2),point(3)) ];
                d(3:4) = [ fun_test_A2(point(1),point(2),point(3)), fun_test_A2p(point(1),point(2),point(3)) ];
                d(5:6) = [ fun_test_A3(point(1),point(2),point(3)), fun_test_A3p(point(1),point(2),point(3)) ];
%                 if d(1)*d(2) >= 1e-3 || d(3)*d(4) >= 1e-3 || d(5)*d(6) >= 1e-3
%                     idx_delete = [idx_delete,i];
%                 end
                if d(1)*d(2) >= 1e-3 || d(3)*d(4) >= 1e-3 || d(5)*d(6) >= 1e-3
                    idx_delete = [idx_delete,i];
                end
            end
            sphere_centers(idx_delete,:) = [];
            sphere_radius(idx_delete)    = [];
            
            sphere_centers       = [sphere_centers;periodic_sphere_centers];
            sphere_radius        = [sphere_radius,periodic_sphere_radius];  
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
            
            n_sphere = length(sphere_radius);
            periodic_sphere_centers       = [];
            periodic_sphere_radius        = [];
            
            for i = 1:size(shift,1)   
                for j = 1:n_sphere
                    center_temp = sphere_centers(j,:) + shift(i,:);
                    point       = [a1,a2,a3]*center_temp';
                    
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
                    vec_n = cross(a1_orig+a2_orig,a3_orig);
                    fun_test_AA13  = @(x,y,z) vec_n(1)*(x-a1_orig(1)) + vec_n(2)*(y-a1_orig(2)) + vec_n(3)*(z-a1_orig(3));
                    fun_test_AA13p = @(x,y,z) vec_n(1)*(x-a2_orig(1)) + vec_n(2)*(y-a2_orig(2)) + vec_n(3)*(z-a2_orig(3));
                    d(5:6) = [ fun_test_AA13(point(1),point(2),point(3)), fun_test_AA13p(point(1),point(2),point(3)) ];
                    if d(1)*d(2) < 0 && d(3)*d(4) < 0 && d(5)*d(6) < 0
                        periodic_sphere_centers = [periodic_sphere_centers; center_temp];
                        periodic_sphere_radius  = [periodic_sphere_radius, sphere_radius(j)];
                    end
                end
            end
            sphere_centers       = [sphere_centers;periodic_sphere_centers];
            sphere_radius        = [sphere_radius,periodic_sphere_radius];
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
            n_sphere = length(sphere_radius);
            periodic_sphere_centers       = [];
            periodic_sphere_radius        = [];
            for i = 1:size(shift,1)   
                for j = 1:n_sphere
                    center_temp = sphere_centers(j,:) + shift(i,:);
                    point       = [a1,a2,a3]*center_temp';
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
                        periodic_sphere_centers = [periodic_sphere_centers; center_temp];
                        periodic_sphere_radius  = [periodic_sphere_radius, sphere_radius(j)];
                    end
                end
            end
            
            idx_delete   = [];
            for i = 1:n_sphere
                center_temp = sphere_centers(i,:);
                point       = [a1,a2,a3]*center_temp';
                d(1:2) = [ fun_test_A1(point(1),point(2),point(3)), fun_test_A1p(point(1),point(2),point(3)) ];
                d(3:4) = [ fun_test_A2(point(1),point(2),point(3)), fun_test_A2p(point(1),point(2),point(3)) ];
                d(5:6) = [ fun_test_A3(point(1),point(2),point(3)), fun_test_A3p(point(1),point(2),point(3)) ];
                if d(1) <= -1e-3 || d(2) >= 1e-3 || d(3) >= 1e-3 || d(4) <= -1e-3 || d(5) <= -1e-3 || d(6) >= 1e-3
                    idx_delete = [idx_delete,i];
                end
            end
            sphere_centers(idx_delete,:) = [];
            sphere_radius(idx_delete)    = [];
            
            sphere_centers       = [sphere_centers;periodic_sphere_centers];
            sphere_radius        = [sphere_radius,periodic_sphere_radius];
        case 'a_base_centered_monoclinic'
            shift = [ 0 1 0;
                      0 0 1;
                      0 1 1];
            n_sphere = length(sphere_radius);
            periodic_sphere_centers       = [];
            periodic_sphere_radius        = [];
            for i = 1:size(shift,1)   
                for j = 1:n_sphere
                    center_temp = sphere_centers(j,:) + shift(i,:);
                    point       = [a1,a2,a3]*center_temp';
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
%                         periodic_sphere_centers = [periodic_sphere_centers; center_temp];
%                         periodic_sphere_radius  = [periodic_sphere_radius, sphere_radius(j)];
%                     end

                    if d(1)*d(2) <= 1e-5 && d(3)*d(4) <= 1e-5 && d(5)*d(6) <= 1e-5
                        periodic_sphere_centers = [periodic_sphere_centers; center_temp];
                        periodic_sphere_radius  = [periodic_sphere_radius, sphere_radius(j)];
                    end
                end
            end
            
            idx_delete   = [];
            for i = 1:n_sphere
                center_temp = sphere_centers(i,:);
                point       = [a1,a2,a3]*center_temp';
                d(1:2) = [ fun_test_A1(point(1),point(2),point(3)), fun_test_A1p(point(1),point(2),point(3)) ];
                d(3:4) = [ fun_test_A2(point(1),point(2),point(3)), fun_test_A2p(point(1),point(2),point(3)) ];
                d(5:6) = [ fun_test_A3(point(1),point(2),point(3)), fun_test_A3p(point(1),point(2),point(3)) ];
%                 if d(1)*d(2) >= 1e-3 || d(3)*d(4) >= 1e-3 || d(5)*d(6) >= 1e-3
%                     idx_delete = [idx_delete,i];
%                 end
                if d(1)*d(2) >= -1e-5 || d(3)*d(4) >= -1e-5 || d(5)*d(6) >= -1e-5
                    idx_delete = [idx_delete,i];
                end
            end
            sphere_centers(idx_delete,:) = [];
            sphere_radius(idx_delete)    = [];
            
            sphere_centers       = [sphere_centers;periodic_sphere_centers];
            sphere_radius        = [sphere_radius,periodic_sphere_radius];    
        case 'a_base_centered_orthorhombic'
            shift = [ 0 -1 0;
                      0  0 1;
                      0 -1 1];
            n_sphere = length(sphere_radius);
            periodic_sphere_centers       = [];
            periodic_sphere_radius        = [];
            for i = 1:size(shift,1)   
                for j = 1:n_sphere
                    center_temp = sphere_centers(j,:) + shift(i,:);
                    point       = [a1,a2,a3]*center_temp';
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
%                     if d(1)*d(2) <= 1e-3 && d(3)*d(4) <= 1e-3 && d(5)*d(6) <= 1e-3
%                         periodic_sphere_centers = [periodic_sphere_centers; center_temp];
%                         periodic_sphere_radius  = [periodic_sphere_radius, sphere_radius(j)];
%                     end
                    if d(1)*d(2) <= 1 && d(3)*d(4) <= 1 && d(5)*d(6) <= 1
                        periodic_sphere_centers = [periodic_sphere_centers; center_temp];
                        periodic_sphere_radius  = [periodic_sphere_radius, sphere_radius(j)];
                    end
                end
            end
            
            idx_delete   = [];
            for i = 1:n_sphere
                center_temp = sphere_centers(i,:);
                point       = [a1,a2,a3]*center_temp';
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
            sphere_centers(idx_delete,:) = [];
            sphere_radius(idx_delete)    = [];
            
            sphere_centers       = [sphere_centers;periodic_sphere_centers];
            sphere_radius        = [sphere_radius,periodic_sphere_radius];
        case 'c_base_centered_orthorhombic'
            shift = [  0  1 0;
                      -1  0 0;
                      -1  1 0];
            n_sphere = length(sphere_radius);
            periodic_sphere_centers       = [];
            periodic_sphere_radius        = [];
            for i = 1:size(shift,1)   
                for j = 1:n_sphere
                    center_temp = sphere_centers(j,:) + shift(i,:);
                    point       = [a1,a2,a3]*center_temp';
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
%                     if d(1)*d(2) <= 1e-3 && d(3)*d(4) <= 1e-3 && d(5)*d(6) <= 1e-3
%                         periodic_sphere_centers = [periodic_sphere_centers; center_temp];
%                         periodic_sphere_radius  = [periodic_sphere_radius, sphere_radius(j)];
%                     end
                    if d(1)*d(2) <= 1e-5 && d(3)*d(4) <= 1e-5 && d(5)*d(6) <= 1e-5
                        periodic_sphere_centers = [periodic_sphere_centers; center_temp];
                        periodic_sphere_radius  = [periodic_sphere_radius, sphere_radius(j)];
                    end
                end
            end
            
            idx_delete   = [];
            for i = 1:n_sphere
                center_temp = sphere_centers(i,:);
                point       = [a1,a2,a3]*center_temp';
                d(1:2) = [ fun_test_A1(point(1),point(2),point(3)), fun_test_A1p(point(1),point(2),point(3)) ];
                d(3:4) = [ fun_test_A2(point(1),point(2),point(3)), fun_test_A2p(point(1),point(2),point(3)) ];
                d(5:6) = [ fun_test_A3(point(1),point(2),point(3)), fun_test_A3p(point(1),point(2),point(3)) ];
%                 if d(1)*d(2) >= 1e-3 || d(3)*d(4) >= 1e-3 || d(5)*d(6) >= 1e-3
%                     idx_delete = [idx_delete,i];
%                 end
                if d(1)*d(2) >= -1e-5 || d(3)*d(4) >= -1e-5 || d(5)*d(6) >= -1e-5
                    idx_delete = [idx_delete,i];
                end
            end
            sphere_centers(idx_delete,:) = [];
            sphere_radius(idx_delete)    = [];
            
            sphere_centers       = [sphere_centers;periodic_sphere_centers];
            sphere_radius        = [sphere_radius,periodic_sphere_radius];    
    end

end