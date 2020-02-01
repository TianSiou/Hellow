function idx = FAME_Material_Locate_Handle_No221_SC_alpha_Po( X, Y, Z, a1, a2, a3, sphere_radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See detail on following website:
%        http://www.aflowlib.org/CrystalDatabase/A_cP1_221_a.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    % Given sphere centers
    sphere_center_Po     = [0; 0; 0];
   
  
    idx_sphere_Po   =  FAME_Material_Locate_Sphere( [X,Y,Z], sphere_radius, sphere_center_Po' );

 
    
 
    % Merge indices
    idx{1}  =  idx_sphere_Po; % indices inner Po
 
end