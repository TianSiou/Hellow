function FAME_Plot_Material_User_Isosurf(grid_nums, a1, a2, a3, lattice_type, lattice_constant, Omega, isofunction, isovalue, hax, plot_mode)
    
    axes(hax)
    hold on
    %% Compute orginal lattice vectors
    lattice_vec_a        = [a1,a2,a3];
    lattice_vec_a_orig_P = Omega'*lattice_vec_a;
    I = eye(3); P = I(:,lattice_constant.Permutation); invP = P';
    invPermutation = [find(invP(:,1)==1), find(invP(:,2)==1), find(invP(:,3)==1)];
    lattice_vec_a_orig   = lattice_vec_a_orig_P(:,invPermutation);

    invAP  = inv(lattice_vec_a_orig_P);
    %% Plot the frame line of primitive cell
    FAME_Plot_Parallelepiped([0,0,0],lattice_vec_a_orig(:,1),lattice_vec_a_orig(:,2),lattice_vec_a_orig(:,3),lattice_constant,'color','computational');
    
    if strcmp(plot_mode,'unit_cell') == 1
        switch lattice_type
            case {'simple_cubic','triclinic','primitive_monoclinic','primitive_tetragonal','primitive_othorhombic'}
                FAME_Plot_Parallelepiped([0;0;0],a1,a2,a3,lattice_constant,'no_color');
            case {'face_centered_cubic','face_centered_orthorhombic'}
                FAME_Plot_Parallelepiped([0;0;0],-a1+a2+a3,a1-a2+a3,a1+a2-a3,lattice_constant,'no_color');
            case {'body_centered_cubic','body_centered_tetragonal','body_centered_orthorhombic'}
                FAME_Plot_Parallelepiped([0;0;0],a2+a3,a1+a3,a1+a2,lattice_constant,'no_color');
            case 'hexagonal'
                FAME_Plot_Parallelepiped(         a2_orig,         a1_orig, a1_orig+a2_orig, a3_orig ,lattice_constant,'no_color');
                FAME_Plot_Parallelepiped( a1_orig+a2_orig, a1_orig+a2_orig,        -a2_orig, a3_orig, lattice_constant,'no_color');
            case 'rhombohedral'
                FAME_Plot_Parallelepiped([0;0;0],a1-a2,a2-a3,a1+a2+a3,lattice_constant,'no_color');
            case 'a_base_centered_monoclinic'
                FAME_Plot_Parallelepiped(a2,a3-a2,a3+a2,a1,lattice_constant,'no_color');
            case 'a_base_centered_orthorhombic'
                FAME_Plot_Parallelepiped([0;0;0],a2+a3,a3-a2,a1,lattice_constant,'no_color');
            case 'c_base_centered_orthorhombic'
                FAME_Plot_Parallelepiped([0;0;0],a2+a1,a2-a1,a3,lattice_constant,'no_color');    
        end
    end

    %% Construct mesh grid
    temp = [[0;0;0], lattice_vec_a(:,1),lattice_vec_a(:,2),lattice_vec_a(:,3),lattice_vec_a(:,1)+lattice_vec_a(:,2),lattice_vec_a(:,1)+lattice_vec_a(:,3),lattice_vec_a(:,2)+lattice_vec_a(:,3), lattice_vec_a(:,1)+lattice_vec_a(:,2)+lattice_vec_a(:,3)];
    e    = 2/n;
    x_start = min(temp(1,:))-e;  x_end   = max(temp(1,:))+e;
    y_start = min(temp(2,:))-e;  y_end   = max(temp(2,:))+e;
    z_start = min(temp(3,:))-e;  z_end   = max(temp(3,:))+e;
    X = linspace(x_start,x_end,grid_nums(1));  Y = linspace(y_start,y_end,grid_nums(2));  Z = linspace(z_start,z_end,grid_nums(3));
    [X,Y,Z] = meshgrid(X,Y,Z);
    %% Compute inner cell indices
    point_set       = [X(:),Y(:),Z(:)];
    point_set_orig  = point_set*Omega;
    X_orig = reshape( point_set_orig(:,1), grid_nums(1),grid_nums(2),grid_nums(3));  
    Y_orig = reshape( point_set_orig(:,2), grid_nums(1),grid_nums(2),grid_nums(3));  
    Z_orig = reshape( point_set_orig(:,3), grid_nums(1),grid_nums(2),grid_nums(3));
    coef            = point_set_orig*invAP';

    idx_outofcell = union( find(coef(:,1)<0 | coef(:,1)>1),...
                           union( find(coef(:,2)<0 | coef(:,2)>1),find(coef(:,3)<0 | coef(:,3)>1) ) );
    
    %% Plot isosurface    
    color_map = { [15 131 225]/255, [200,40,45]/255, [238,126,2]/255, [75,165,102]/255, rand(1,3), rand(1,3), rand(1,3), rand(1,3), rand(1,3), rand(1,3)};
    G = isofunction(point_set_orig(:,1),point_set_orig(:,2),point_set_orig(:,3),lattice_vec_a_orig(:,1),lattice_vec_a_orig(:,2),lattice_vec_a_orig(:,3));
    for i = 1:length(G)
        G{i}(idx_outofcell) = 0;
        
        p{i} = patch( isosurface(X_orig,Y_orig,Z_orig,reshape(G{i},n,n,n),isovalue(i)) );
        p{i}.FaceColor = color_map{i};
        p{i}.EdgeColor = 'none';
    end
    camlight right
    rotate3d on
    axis equal
    axis off
    view([-38 17])

 end
    