function idx = FAME_Material_Locate_Handle_No216_FCC_AuBe5( X, Y, Z, a1, a2, a3, sphere_radius )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The lattice constant for this structure is taken from (Batchelder, 1958),which does not give the internal coordinate for the (16c) site.
%However, (Baenziger, 1950) assumes that uranium compounds of this type have an internal parameter x3?5/8.(Pearson, 1958) uses this to infer a value of x3?5/8 here as well.
% See detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/AB5_cF24_216_a_ce.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x3 = 5/8;
    % Given sphere centers
    sphere_center_Au     =  [0; 0; 0];
    sphere_center_Be_1  =    (1/4)*a1 +    (1/4)*a2 +    (1/4)*a3;
    sphere_center_Be_2  =       x3*a1 +       x3*a2 +       x3*a3;
    sphere_center_Be_3  =       x3*a1 +       x3*a2 + (2-3*x3)*a3;
    sphere_center_Be_4  =       x3*a1 + (2-3*x3)*a2 +       x3*a3;
    sphere_center_Be_5  = (2-3*x3)*a1 +       x3*a2 +       x3*a3;
    % Determine indices inner material
    idx_sphere_Be_1  =  FAME_Material_Locate_Sphere( [X,Y,Z], sphere_radius, sphere_center_Be_1' );
    idx_sphere_Be_2  =  FAME_Material_Locate_Sphere( [X,Y,Z], sphere_radius, sphere_center_Be_2' );
    idx_sphere_Be_3  =  FAME_Material_Locate_Sphere( [X,Y,Z], sphere_radius, sphere_center_Be_3' );
    idx_sphere_Be_4  =  FAME_Material_Locate_Sphere( [X,Y,Z], sphere_radius, sphere_center_Be_4' );
    idx_sphere_Be_5  =  FAME_Material_Locate_Sphere( [X,Y,Z], sphere_radius, sphere_center_Be_5' );
   
    
    idx_sphere_Au   =  FAME_Material_Locate_Sphere( [X,Y,Z], sphere_radius, sphere_center_Au' );
    idx_sphere_Be  =  union(idx_sphere_Be_1,union(idx_sphere_Be_2,union(idx_sphere_Be_3,union(idx_sphere_Be_4,idx_sphere_Be_5)))) ;
    
   
    
    % Merge indices
    idx{1}  = idx_sphere_Au; % indices inner Au
    idx{2}  = idx_sphere_Be; % indices inner Be
end