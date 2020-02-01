function [ Par_mesh, Par_lattice, Par_recip_lattice, Par_material, Par_eig ] = FAME_Parameter_Generator( Popt )
% This routine return the computational informations about given input:
%   grid_nums: a 1x3 array which contains grid numbers in terms of x-,y-, and z-axes, these must be positive integers.
%   lattice type: 
%   lattice_constant:
%   Par_material:

    grid_num     = Popt.mesh.grid_num;
    lattice_type = Popt.lattice.lattice_type; 
    Par_material = Popt.material;
    lattice_constant = [];
    if isfield( Popt.lattice, 'lattice_constant')
        lattice_constant_orig = Popt.lattice.lattice_constant; 
    end
    if isfield( Popt.lattice, 'lattice_vec_a')
        lattice_vec_a_orig = Popt.lattice.lattice_vec_a; 
    end

    if strcmp(lattice_type,'user_defined')
        % read material infomations from 'Par_material.data_name'
        data         = feval(Par_material.data_name,Par_material.sphere_radius,Par_material.cylinder_radius);
        lattice_type = data.lattice_type;
        Par_material.material_handle   = @(x,y,z,a1,a2,a3) FAME_Material_Locate_Handle_User_Defined( x, y, z, a1, a2, a3, Par_material.data_name,Par_material.sphere_radius,Par_material.cylinder_radius);
        
        if isstruct(lattice_constant) == 0
            lattice_constant_orig = data.lattice_constant;
        end
        Par_material.data = data;
    elseif isfield(Popt.material, 'isofunction')
        Par_material.material_handle   = @(x,y,z,a1,a2,a3) FAME_Material_Locate_Handle_User_Defined( x, y, z, a1, a2, a3, Popt.material.isofunction, Popt.material.isovalue);
    end
    
            
    switch lattice_type
        case {'simple_cubic', 'primitive_orthorhombic', 'primitive_tetragonal'}
            lattice_constant_orig = FAME_Parameter_Lattice_Constants_Format(lattice_type,lattice_constant_orig);
            
            lattice_constant.Permutation = [1 2 3];
            Omega = eye(3);

            lattice_vec_a1 = lattice_constant_orig.a * [1;0;0];
            lattice_vec_a2 = lattice_constant_orig.b * [0;1;0];
            lattice_vec_a3 = lattice_constant_orig.c * [0;0;1];
            lattice_vec_a  = [lattice_vec_a1,lattice_vec_a2,lattice_vec_a3];
            lattice_vec_a_orig = lattice_vec_a;
            
            edge_len(1) = lattice_constant_orig.a;
            edge_len(2) = lattice_constant_orig.b;
            edge_len(3) = lattice_constant_orig.c;

            mesh_len     = [edge_len(1)/grid_num(1) edge_len(2)/grid_num(2) edge_len(3)/grid_num(3)];
        otherwise
            if isfield( Popt.lattice, 'lattice_vec_a') ~= 1
                lattice_constant_orig = FAME_Parameter_Lattice_Constants_Format(lattice_type,lattice_constant_orig);
                
                [lattice_vec_a,lattice_vec_a_orig,lattice_constant.length_a1,lattice_constant.length_a2,lattice_constant.length_a3,lattice_constant.theta_1,lattice_constant.theta_2,lattice_constant.theta_3,lattice_constant.Permutation] =...
                    FAME_Parameter_Lattice_Vector(lattice_type,lattice_constant_orig);
            else
                lattice_constant_orig = FAME_Parameter_Lattice_Constants_Format(lattice_type,lattice_vec_a_orig);
                
                [lattice_vec_a,lattice_vec_a_orig,lattice_constant.length_a1,lattice_constant.length_a2,lattice_constant.length_a3,lattice_constant.theta_1,lattice_constant.theta_2,lattice_constant.theta_3,lattice_constant.Permutation] =...
                    FAME_Parameter_Lattice_Vector(lattice_type,lattice_constant_orig,lattice_vec_a_orig);
            end
            
            Omega = lattice_vec_a/lattice_vec_a_orig(:,lattice_constant.Permutation);
            I = eye(3); P = I(:,lattice_constant.Permutation); invP = P';
            lattice_constant.invPermutation = [find(invP(:,1)==1), find(invP(:,2)==1), find(invP(:,3)==1)];
            
            edge_len(1) = lattice_vec_a(1,1);
            edge_len(2) = lattice_vec_a(2,2);
            edge_len(3) = lattice_vec_a(3,3);              
            mesh_len    = [edge_len(1)/grid_num(1) edge_len(2)/grid_num(2) edge_len(3)/grid_num(3)];

            [lattice_constant.m1,lattice_constant.m2,lattice_constant.m3,lattice_constant.m4,lattice_constant.t1,lattice_constant.t2,lattice_constant.t3,lattice_constant.t4] = ...
                FAME_Parameter_Boundary_Point(lattice_vec_a,lattice_constant,mesh_len,grid_num);
    end
    
    Par_mesh.grid_num             = grid_num;
    Par_mesh.edge_len             = edge_len;
    Par_mesh.mesh_len             = mesh_len;
    Par_lattice.lattice_type      = lattice_type;
    Par_lattice.lattice_constant  = lattice_constant;
    Par_lattice.lattice_constant_orig  = lattice_constant_orig;
    Par_lattice.lattice_vec_a     = lattice_vec_a;
    Par_lattice.lattice_vec_a_orig = lattice_vec_a_orig;
    Par_recip_lattice.reciprocal_lattice_vector_b = inv(lattice_vec_a');
    Par_recip_lattice.reciprocal_lattice_vector_b_orig = inv(lattice_vec_a_orig');
    Par_lattice.Omega             = Omega;
    Par_eig                       = Popt.eig;
end


% function data = periodic_check_simple( data, lattice_vec_a )
%     
%     a = lattice_vec_a(1,1); b = lattice_vec_a(2,2); c = lattice_vec_a(3,3);
%     for i = 1:length(data.parameters)
%         sphere_centers = data.parameters{i}.sphere_centers;
%         sphere_radius  = data.parameters{i}.sphere_radius(1);
%         cylinder_bot_centers = data.parameters{i}.cylinder_bot_centers;
%         cylinder_top_centers = data.parameters{i}.cylinder_top_centers;
%         cylinder_radius      = data.parameters{i}.cylinder_radius(1);
%         
%         sphere_centers = [0,0,0;0.5,0,0;1,0,0];
%         %% For sphere
%         % 週期生成 27 組 sphere_centers
%         e = ones(3,1);  s = [-1;0;1];
%         shift = [kron(kron(e,e),s), kron(kron(e,s),e), kron(kron(s,e),e)];
%         e = ones(size(sphere_centers,1),1);
%         shift_tmp = kron(e,shift);
%         
%         e = ones(size(shift,1),1);
%         sphere_centers_tmp = kron(sphere_centers,e);
%         sphere_centers_tmp = sphere_centers_tmp + shift_tmp;
%         % 排除相同的 sphere_center
%         tol = 1e-6;
%         for j = 1:size(sphere_centers_tmp,1)-1
%             center = sphere_centers_tmp(j,:);
%             idx_same = find( abs(sphere_centers_tmp(j+1:end,1) - center(1)) < tol & ...
%                              abs(sphere_centers_tmp(j+1:end,2) - center(2)) < tol & ...
%                              abs(sphere_centers_tmp(j+1:end,3) - center(3)) < tol );
%             idx_same = idx_same + j;
%             sphere_centers_tmp(idx_same,:) = nan;
%         end
%         idx_nan = find(isnan(sphere_centers_tmp(:,1)));
%         sphere_centers_tmp(idx_nan,:) = [];
%         % 判斷球與 primitive cell 是否有相交
%         idx_outofcell = union(union( find( sphere_centers_tmp(:,1) > 1 | sphere_centers_tmp(:,1) < 0 ),...
%                                      find( sphere_centers_tmp(:,2) > 1 | sphere_centers_tmp(:,2) < 0 ) ),...
%                                      find( sphere_centers_tmp(:,3) > 1 | sphere_centers_tmp(:,3) < 0 )   );
%         idx_delete = [];
%         for j = 1:length(idx_outofcell)
%             center       = sphere_centers_tmp(idx_outofcell(j),:);
%             exact_center = center*lattice_vec_a';
%             
%             flag_a1 = ( center(1) > 1 && (exact_center(1) - a > sphere_radius) ) || ( center(1) < 0 && (0 - exact_center(1) > sphere_radius) );
%             flag_a2 = ( center(2) > 1 && (exact_center(2) - b > sphere_radius) ) || ( center(2) < 0 && (0 - exact_center(2) > sphere_radius) );
%             flag_a3 = ( center(3) > 1 && (exact_center(3) - c > sphere_radius) ) || ( center(3) < 0 && (0 - exact_center(3) > sphere_radius) );
% 
%             if (flag_a1 || flag_a2 || flag_a3) == 1
%                 idx_delete = [idx_delete, idx_outofcell(j)];
%             end
%             
%         end
%         sphere_centers_tmp(idx_delete,:) = [];
%         data.parameters{i}.sphere_centers = sphere_centers_tmp;
%         %% For cylinder
%         % 週期生成 27 組 cylinder_bot_centers 與 cylinder_top_centers
%         e = ones(3,1);  s = [-1;0;1];
%         shift = [kron(kron(e,e),s), kron(kron(e,s),e), kron(kron(s,e),e)];
%         e = ones(size(cylinder_bot_centers,1),1);
%         shift_tmp = kron(e,shift);
%         
%         e = ones(size(shift,1),1);
%         cylinder_bot_centers_tmp = kron(cylinder_bot_centers,e);
%         cylinder_bot_centers_tmp = cylinder_bot_centers_tmp + shift_tmp;
%         cylinder_top_centers_tmp = kron(cylinder_top_centers,e);
%         cylinder_top_centers_tmp = cylinder_top_centers_tmp + shift_tmp;
%         % 排除相同的 cylinder_bot_centers_tmp, cylinder_top_centers_tmp 組
%         tol = 1e-6;
%         for j = 1:size(cylinder_bot_centers_tmp,1)-1
%             bot_center = cylinder_bot_centers_tmp(j,:);
%             top_center = cylinder_top_centers_tmp(j,:);
%             idx_same_bot = find( abs(cylinder_bot_centers_tmp(j+1:end,1) - bot_center(1)) < tol & ...
%                                  abs(cylinder_bot_centers_tmp(j+1:end,2) - bot_center(2)) < tol & ...
%                                  abs(cylinder_bot_centers_tmp(j+1:end,3) - bot_center(3)) < tol );
%             idx_same_top = find( abs(cylinder_top_centers_tmp(j+1:end,1) - top_center(1)) < tol & ...
%                                  abs(cylinder_top_centers_tmp(j+1:end,2) - top_center(2)) < tol & ...
%                                  abs(cylinder_top_centers_tmp(j+1:end,3) - top_center(3)) < tol );
%             idx_same_bot = idx_same_bot + j;
%             idx_same_top = idx_same_top + j;
%             idx_same     = intersect(idx_same_bot,idx_same_top);
%             cylinder_bot_centers_tmp(idx_same,:) = nan;
%             cylinder_top_centers_tmp(idx_same,:) = nan;
%         end
%         idx_nan = find(isnan(cylinder_bot_centers_tmp(:,1)));
%         cylinder_bot_centers_tmp(idx_nan,:) = [];
%         cylinder_top_centers_tmp(idx_nan,:) = [];
%         % 判斷圓柱與 primitive cell 六個面是否有相交
%         idx_innerofcell = [];
%         for j = 1:size(cylinder_bot_centers_tmp,1)
%             bot_center = cylinder_bot_centers_tmp(j,:);
%             top_center = cylinder_top_centers_tmp(j,:);
%             exact_bot_center = bot_center*lattice_vec_a';
%             exact_top_center = top_center*lattice_vec_a';
%             cylinder_vec     = exact_top_center-exact_bot_center;
%             
%             
%             if abs( dot(cylinder_vec',cross(lattice_vec_a(:,1), lattice_vec_a(:,2))) ) > tol
%                 % 與 (原點-a1-a2) 面 / (a3點-a1-a2) 面
%                 A_o_12  = [ cylinder_vec', lattice_vec_a(:,1), lattice_vec_a(:,2) ];
%                 A_a3_12 = [ cylinder_vec', lattice_vec_a(:,1)-lattice_vec_a(:,3), lattice_vec_a(:,2)-lattice_vec_a(:,3) ];
%                 tmp_o_12  = A_o_12\top_center';
%                 tmp_a3_12 = A_a3_12\(top_center'-lattice_vec_a(:,3));
%                 t_o_12  = tmp_o_12(1);  u_o_12  = tmp_o_12(2);  v_o_12  = tmp_o_12(3);
%                 t_a3_12 = tmp_a3_12(1); u_a3_12 = tmp_a3_12(2); v_a3_12 = tmp_a3_12(3);
%             else
%                 t_o_12  = nan; u_o_12  = nan; v_o_12  = nan;
%                 t_a3_12 = nan; u_a3_12 = nan; v_a3_12 = nan;
%             end
%             if abs( dot(cylinder_vec',cross(lattice_vec_a(:,1), lattice_vec_a(:,3))) ) > tol
%                 % 與 (原點-a1-a3) 面 / (a2點-a1-a3) 面
%                 A_o_13  = [ cylinder_vec', lattice_vec_a(:,1), lattice_vec_a(:,3) ];
%                 A_a2_13 = [ cylinder_vec', lattice_vec_a(:,1)-lattice_vec_a(:,2), lattice_vec_a(:,3)-lattice_vec_a(:,2) ];
%                 tmp_o_13  = A_o_13\top_center';
%                 tmp_a2_13 = A_a2_13\(top_center'-lattice_vec_a(:,2));
%                 t_o_13  = tmp_o_13(1);  u_o_13  = tmp_o_13(2);  v_o_13  = tmp_o_13(3);
%                 t_a2_13 = tmp_a2_13(1); u_a2_13 = tmp_a2_13(2); v_a2_13 = tmp_a2_13(3);
%             else
%                 t_o_13  = nan; u_o_13  = nan; v_o_13  = nan;
%                 t_a2_13 = nan; u_a2_13 = nan; v_a2_13 = nan;
%             end
%             if abs( dot(cylinder_vec',cross(lattice_vec_a(:,2), lattice_vec_a(:,3))) ) > tol
%                 % 與 (原點-a2-a3) 面 / (a1點-a2-a3) 面
%                 A_o_23  = [ cylinder_vec', lattice_vec_a(:,2), lattice_vec_a(:,3) ];
%                 A_a1_23 = [ cylinder_vec', lattice_vec_a(:,2)-lattice_vec_a(:,1), lattice_vec_a(:,3)-lattice_vec_a(:,1) ];
%                 tmp_o_23  = A_o_23\top_center';
%                 tmp_a1_23 = A_a1_23\(top_center'-lattice_vec_a(:,1));
%                 t_o_23  = tmp_o_23(1);  u_o_23  = tmp_o_23(2);  v_o_23  = tmp_o_23(3);
%                 t_a1_23 = tmp_a1_23(1); u_a1_23 = tmp_a1_23(2); v_a1_23 = tmp_a1_23(3);
%             else
%                 t_o_23  = nan; u_o_23  = nan; v_o_23  = nan;
%                 t_a1_23 = nan; u_a1_23 = nan; v_a1_23 = nan;
%             end
%             
%             flag_o_12  = (t_o_12 >= 0 && t_o_12 < 1 ) && (u_o_12 >= 0 && u_o_12 <= 1 ) && (v_o_12 >= 0 && v_o_12 <= 1 ); % 表示與 (原點-a1-a2) 面 相交
%             flag_a3_12 = (t_a3_12 >= 0 && t_a3_12 < 1 ) && (u_a3_12 >= 0 && u_a3_12 <= 1 ) && (v_a3_12 >= 0 && v_a3_12 <= 1 ); % 表示與 (a3點-a1-a2) 面 相交
%             flag_o_13  = (t_o_13 >= 0 && t_o_13 < 1 ) && (u_o_13 >= 0 && u_o_13 <= 1 ) && (v_o_13 >= 0 && v_o_13 <= 1 ); % 表示與 (原點-a1-a3) 面 相交
%             flag_a2_13 = (t_a2_13 >= 0 && t_a2_13 < 1 ) && (u_a2_13 >= 0 && u_a2_13 <= 1 ) && (v_a2_13 >= 0 && v_a2_13 <= 1 ); % 表示與 (a2點-a1-a3) 面 相交
%             flag_o_23  = (t_o_23 >= 0 && t_o_23 < 1 ) && (u_o_23 >= 0 && u_o_23 <= 1 ) && (v_o_23 >= 0 && v_o_23 <= 1 ); % 表示與 (原點-a2-a3) 面 相交
%             flag_a1_23 = (t_a1_23 >= 0 && t_a1_23 < 1 ) && (u_a1_23 >= 0 && u_a1_23 <= 1 ) && (v_a1_23 >= 0 && v_a1_23 <= 1 ); % 表示與 (a1點-a2-a3) 面 相交
%             flag_inner = flag_o_12 || flag_a3_12 || flag_o_13 || flag_a2_13 || flag_o_23 || flag_a1_23 ;
%             if flag_inner == 1
%                 idx_innerofcell = [idx_innerofcell, j];
%             end
%         end
%         test1 = cylinder_bot_centers_tmp(idx_innerofcell,:);
%         test2 = cylinder_top_centers_tmp(idx_innerofcell,:);
%     end
% end
