function idx = FAME_Material_Locate_Handle_No198_SC_alpha_CO( X, Y, Z, a1, a2, a3, sphere_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/AB_cP8_198_a_a.alpha-CO.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
     % Given sphere centers
    x1 = 0.15;
    x2=  0.85;
    sphere_center_C_1  =       x1*a1  +       x1*a2 +        x1*a3;
    sphere_center_C_2  = (1/2-x1)*a1  +   (1-x1)*a2 +  (1/2+x1)*a3;
    sphere_center_C_3  =  (1 -x1)*a1  + (1/2+x1)*a2 +  (1/2-x1)*a3;
    sphere_center_C_4  = (1/2+x1)*a1  + (1/2-x1)*a2 +    (1-x1)*a3;
    sphere_center_O_1  =       x2*a1  +       x2*a2 +        x2*a3;
    sphere_center_O_2  = (1/2-x2)*a1  +   (1-x2)*a2 +  (1/2+x2)*a3;
    sphere_center_O_3  =   (1-x2)*a1  + (1/2+x2)*a2 +  (1/2-x2)*a3;
    sphere_center_O_4  = (1/2+x2)*a1  + (1/2-x2)*a2 +    (1-x2)*a3;
    
    % Determine indices inner material
    idx_sphere_C_1  =  FAME_Material_Locate_Sphere( [X,Y,Z], sphere_radius, sphere_center_O_1' );
    idx_sphere_C_2  =  FAME_Material_Locate_Sphere( [X,Y,Z], sphere_radius, sphere_center_O_2' );
    idx_sphere_C_3  =  FAME_Material_Locate_Sphere( [X,Y,Z], sphere_radius, sphere_center_O_3' );
    idx_sphere_C_4  =  FAME_Material_Locate_Sphere( [X,Y,Z], sphere_radius, sphere_center_O_4' );
    idx_sphere_O_1  =  FAME_Material_Locate_Sphere( [X,Y,Z], sphere_radius, sphere_center_C_1' );
    idx_sphere_O_2  =  FAME_Material_Locate_Sphere( [X,Y,Z], sphere_radius, sphere_center_C_2' );
    idx_sphere_O_3  =  FAME_Material_Locate_Sphere( [X,Y,Z], sphere_radius, sphere_center_C_3' );
    idx_sphere_O_4  =  FAME_Material_Locate_Sphere( [X,Y,Z], sphere_radius, sphere_center_C_4' );
    
   
    idx_sphere_C  =  union(idx_sphere_C_1,union(idx_sphere_C_2,union(idx_sphere_C_3,idx_sphere_C_4)));
    idx_sphere_O  =  union(idx_sphere_O_1,union(idx_sphere_O_2,union(idx_sphere_O_3,idx_sphere_O_4)));
    
 
    % Merge indices
    idx{1}  =  idx_sphere_C; % indices inner C
    idx{2}  =  idx_sphere_O; % indices inner O
 
end