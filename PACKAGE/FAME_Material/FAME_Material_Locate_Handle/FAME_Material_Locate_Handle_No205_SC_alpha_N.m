function idx = FAME_Material_Locate_Handle_No205_SC_alpha_N( X, Y, Z, a1, a2, a3, sphere_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/A_cP8_205_c.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
     % Given sphere centers
    x1 = 0.15;
    sphere_center_N_1  =       x1*a1  +       x1*a2 +        x1*a3;
    sphere_center_N_2  = (1/2-x1)*a1  +   (1-x1)*a2 +  (1/2+x1)*a3;
    sphere_center_N_3  =  (1 -x1)*a1  + (1/2+x1)*a2 +  (1/2-x1)*a3;
    sphere_center_N_4  = (1/2+x1)*a1  + (1/2-x1)*a2 +    (1-x1)*a3;
    sphere_center_N_5  =   (1-x1)*a1  +   (1-x1)*a2 +    (1-x1)*a3;
    sphere_center_N_6  = (1/2+x1)*a1  +       x1*a2 +  (1/2-x1)*a3;
    sphere_center_N_7  =       x1*a1  + (1/2-x1)*a2 +  (1/2+x1)*a3;
    sphere_center_N_8  = (1/2-x1)*a1  + (1/2+x1)*a2 +        x1*a3;
    
    % Determine indices inner material
    idx_sphere_N_1  =  FAME_Material_Locate_Sphere( [X,Y,Z], sphere_radius, sphere_center_N_1' );
    idx_sphere_N_2  =  FAME_Material_Locate_Sphere( [X,Y,Z], sphere_radius, sphere_center_N_2' );
    idx_sphere_N_3  =  FAME_Material_Locate_Sphere( [X,Y,Z], sphere_radius, sphere_center_N_3' );
    idx_sphere_N_4  =  FAME_Material_Locate_Sphere( [X,Y,Z], sphere_radius, sphere_center_N_4' );
    idx_sphere_N_5  =  FAME_Material_Locate_Sphere( [X,Y,Z], sphere_radius, sphere_center_N_5' );
    idx_sphere_N_6  =  FAME_Material_Locate_Sphere( [X,Y,Z], sphere_radius, sphere_center_N_6' );
    idx_sphere_N_7  =  FAME_Material_Locate_Sphere( [X,Y,Z], sphere_radius, sphere_center_N_7' );
    idx_sphere_N_8  =  FAME_Material_Locate_Sphere( [X,Y,Z], sphere_radius, sphere_center_N_8' );
    
   
    idx_sphere_N  =  union(idx_sphere_N_1,union(idx_sphere_N_2,union(idx_sphere_N_3,union(idx_sphere_N_4,union(idx_sphere_N_5,union(idx_sphere_N_6,union(idx_sphere_N_7,idx_sphere_N_8)))))));
 
    
 
    % Merge indices
    idx{1}  =  idx_sphere_N; % indices inner N
 
end